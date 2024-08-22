import array
import gzip
import json
import logging
import os
import sqlite3 as sqlite
import subprocess
from math import lgamma
from tempfile import TemporaryDirectory
from typing import Dict, Generator, List, Optional, Tuple

import numpy as np
import scipy.sparse as sparse
from numba import jit

from .cell_calling import call_cells
from .loompy import connect, create


# Copied from cytograph
@jit("float32(float64[:], float64[:])", nopython=True, parallel=True, nogil=True)
def multinomial_distance(p: np.ndarray, q: np.ndarray) -> float:
	N = p.shape[0]
	p_sum = p.sum()
	q_sum = q.sum()
	x = lgamma(N) + lgamma(p_sum + q_sum + N) - lgamma(p_sum + N) - lgamma(q_sum + N)
	for k in range(N):
		x += lgamma(p[k] + 1) + lgamma(q[k] + 1) - lgamma(1) - lgamma(p[k] + q[k] + 1)
	x = np.exp(x)
	return 1 - 1 / (1 + x)


# https://maciejkula.github.io/2015/02/22/incremental-construction-of-sparse-matrices/
class IncrementalSparseMatrixUInt16:
	def __init__(self, shape: Tuple[int, int]):
		self.dtype = np.uint16
		self.shape = shape

		self.rows = array.array('i')
		self.cols = array.array('i')
		self.data = array.array('H')

	def append(self, i: int, j: int, v: int) -> None:
		m, n = self.shape

		if (i >= m or j >= n):
			raise Exception('Index out of bounds')

		self.rows.append(i)
		self.cols.append(j)
		self.data.append(v)

	def tocoo(self) -> sparse.coo_matrix:
		rows = np.frombuffer(self.rows, dtype=np.int32)
		cols = np.frombuffer(self.cols, dtype=np.int32)
		data = np.frombuffer(self.data, dtype=np.uint16)
		return sparse.coo_matrix((data, (rows, cols)), shape=self.shape)

	def __len__(self) -> int:
		return len(self.data)


twobit_to_dna_table = {0: "A", 1: "C", 2: "G", 3: "T"}
dna_to_twobit_table = {"A": 0, "C": 1, "G": 2, "T": 3}


@jit(nopython=True)
def twobit_to_dna(twobit: int, size: int) -> str:
	result = []
	for i in range(size):
		x = (twobit & (3 << 2 * i)) >> 2 * i
		if x == 0:
			result.append("A")
		elif x == 1:
			result.append("C")
		elif x == 2:
			result.append("G")
		elif x == 3:
			result.append("T")
	result.reverse()
	return "".join(result)


@jit(nopython=True)
def dna_to_twobit(dna: str) -> int:
	x = 0
	for nt in dna:
		if nt == "A":
			x += 0
		elif nt == "C":
			x += 1
		elif nt == "G":
			x += 2
		elif nt == "T":
			x += 3
		x <<= 2
	x >>= 2
	return x


@jit(nopython=True)
def twobit_1hamming(twobit: int, size: int) -> List[int]:
	result = []
	for i in range(size):
		x = (twobit >> 2 * (size - i - 1)) & 3
		for j in range(4):
			if x == j:
				continue
			result.append(twobit & ~(3 << 2 * (size - i - 1)) | (j << 2 * (size - i - 1)))
	return result


def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool = True) -> np.ndarray:
	"This is super duper magic sauce to make the order of one list to be like another"
	if check_content:
		assert len(np.intersect1d(a, b)) == len(a), f"The two arrays are not matching"
	return np.argsort(a)[np.argsort(np.argsort(b))]


# TODO: This function is a copy of the same function in loompy.metadata_loaders, call that one instead

def load_sample_metadata(path: str, sample_id: str) -> Dict[str, str]:
	if not os.path.exists(path):
		raise ValueError(f"Samples metadata file '{path}' not found.")
	if path.endswith(".db"):
		# sqlite3
		with sqlite.connect(path) as db:
			cursor = db.cursor()
			cursor.execute("SELECT * FROM sample WHERE name = ?", (sample_id,))
			keys = [x[0] for x in cursor.description]
			vals = cursor.fetchone()
			if vals is not None:
				return dict(zip(keys, vals))
			raise ValueError(f"SampleID '{sample_id}' was not found in the samples database.")
	else:
		result = {}
		with open(path) as f:
			headers = [x.lower() for x in f.readline()[:-1].split("\t")]
			if "sampleid" not in headers and 'name' not in headers:
				raise ValueError("Required column 'SampleID' or 'Name' not found in sample metadata file")
			if "sampleid" in headers:
				sample_metadata_key_idx = headers.index("sampleid")
			else:
				sample_metadata_key_idx = headers.index("name")
			sample_found = False
			for line in f:
				items = line[:-1].split("\t")
				if len(items) > sample_metadata_key_idx and items[sample_metadata_key_idx] == sample_id:
					for i, item in enumerate(items):
						result[headers[i]] = item
					sample_found = True
		if not sample_found:
			raise ValueError(f"SampleID '{sample_id}' not found in sample metadata file")
		return result


class BusFile:
	def __init__(self, path: str, genes_metadata_file: str, genes_metadata_key: str, fragments2genes_file: str, equivalence_classes_file: str, fragments_file: str) -> None:
		self.matrix: sparse.coo_matrix = None
		logging.info("Loading gene metadata")
		self.genes: Dict[str, List[str]] = {}  # Keys are Accessions, values are lists of attribute values
		self.gene_metadata_attributes: List[str] = []  # Attribute names
		with open(genes_metadata_file) as f:
			line = f.readline()
			self.gene_metadata_attributes = line[:-1].split("\t")
			if genes_metadata_key not in self.gene_metadata_attributes:
				raise ValueError(f"Metadata key '{genes_metadata_key}' not found in gene metadata file")
			key_col = self.gene_metadata_attributes.index(genes_metadata_key)
			for line in f:
				items = line[:-1].split("\t")
				self.genes[items[key_col]] = items
		self.accessions = np.array([x for x in self.genes.keys()])
		accession_idx = {acc: i for (i, acc) in enumerate(self.accessions)}
		self.n_genes = len(self.accessions)

		logging.info("Loading fragments-to-gene mappings")
		self.gene_for_fragment: List[str] = []
		with open(fragments2genes_file) as f:
			for line in f:
				transcript_id, accession = line[:-1].split("\t")
				self.gene_for_fragment.append(accession)
		self.gene_for_fragment = np.array(self.gene_for_fragment)

		logging.info("Indexing genes")
		# Array of indices into self.accessions for each gene in gene_for_transcript
		self.gene_for_fragment_idx = np.zeros(len(self.gene_for_fragment), dtype="int32")
		for i in range(len(self.gene_for_fragment)):
			self.gene_for_fragment_idx[i] = accession_idx[self.gene_for_fragment[i]]

		logging.info("Loading equivalence classes")
		self.equivalence_classes: Dict[int, List[int]] = {}
		with open(equivalence_classes_file) as f:
			for line in f:
				ec, trs = line[:-1].split("\t")
				# Each equivalence class is a set of fragments (transcripts)
				self.equivalence_classes[int(ec)] = np.array([int(x) for x in trs.split(",")])
		# But we want each equivalence class mapped to a gene (or -1 if multimapping)
		logging.info("Mapping equivalence classes to genes")
		self.gene_for_ec: Dict[int, int] = {}
		for eqc in self.equivalence_classes.keys():
			gene = -1
			for tid in self.equivalence_classes[eqc]:
				if gene == -1:
					gene = self.gene_for_fragment_idx[tid]
					continue
				elif self.gene_for_fragment_idx[tid] != gene:
					gene = -1  # Multimapping UMI
					break
			self.gene_for_ec[eqc] = gene

		logging.info("Loading fragment IDs")
		self.fragments: List[str] = []
		with open(fragments_file) as f:
			for line in f:
				self.fragments.append(line[:-1])
		self.fragments = np.array(self.fragments)

		logging.info("Loading BUS records")
		fsize = os.path.getsize(path)
		with open(path, "rb") as fb:
			# Read the header
			magic = fb.read(4)
			if magic != b"BUS\0":
				raise IOError("Not a valid BUS file (four leading magic bytes are missing)")
			self.version = int.from_bytes(fb.read(4), byteorder="little", signed=False)
			self.barcode_length = int.from_bytes(fb.read(4), byteorder="little", signed=False)
			self.umi_length = int.from_bytes(fb.read(4), byteorder="little", signed=False)
			tlen = int.from_bytes(fb.read(4), byteorder="little", signed=False)
			self.header = fb.read(tlen).decode("utf-8")  # BUS does not specify an encoding, but let's assume UTF8

			# Read the records
			self.n_records = (fsize - tlen - 20) // 32
			self.bus = np.fromfile(fb, dtype=[
				("barcode", np.uint64),
				("UMI", np.uint64),
				("equivalence_class", np.int32),
				("count", np.uint32),
				("flags", np.uint32),
				("padding", np.int32)
			], count=self.n_records)
			self.bus_gene = np.array([self.gene_for_ec[x] for x in self.bus["equivalence_class"]], dtype=np.int32)
			self.bus_valid = self.bus_gene != -1

		logging.info("Sorting cell IDs")
		self.cell_ids = np.unique(self.bus["barcode"])
		self.cell_ids.sort()
		self.n_cells = len(self.cell_ids)  # This will change after error correction
		self.layers: Dict[str, sparse.coo_matrix] = {}  # Dict of layer name -> sparse matrix
		self.ambient_umis = 0

	def correct(self, whitelist_file: str = None) -> np.ndarray:
		if whitelist_file is not None:
			size = 0
			whitelist = set()
			with open(whitelist_file) as f:
				for bc in f:
					size = len(bc) - 1  # Don't count the newline
					whitelist.add(dna_to_twobit(bc[:-1]))
			for i, bc in enumerate(self.bus["barcode"]):
				if bc in whitelist:
					continue
				corrected = False
				for mut in twobit_1hamming(int(bc), size=size):
					if mut in whitelist:
						self.bus["barcode"][i] = mut
						corrected = True
						break
				if not corrected:
					self.bus_valid[i] = False
		self.cell_ids = np.unique(self.bus["barcode"][self.bus_valid])
		self.cell_ids.sort()
		self.cell_for_barcode_idx = {bc: i for (i, bc) in enumerate(self.cell_ids)}
		self.n_cells = len(self.cell_ids)

	def deduplicate(self) -> None:
		# Sort by barcode, then by UMI, then by gene
		ordering = np.lexsort((self.bus_gene, self.bus["UMI"], self.bus["barcode"]))
		self.bus = self.bus[ordering]
		self.bus_gene = self.bus_gene[ordering]
		self.bus_valid = self.bus_valid[ordering]
		dupes = (self.bus["barcode"][1:] == self.bus["barcode"][:-1]) & (self.bus["UMI"][1:] == self.bus["UMI"][:-1]) & (self.bus_gene[1:] == self.bus_gene[:-1])
		self.bus_valid[1:][dupes] = False

	def count(self) -> sparse.coo_matrix:
		logging.info("Counting pseudoalignments for main matrix")
		genes = self.bus_gene[self.bus_valid]
		cells = [self.cell_for_barcode_idx[x] for x in self.bus["barcode"][self.bus_valid]]
		self.matrix = sparse.coo_matrix((np.ones_like(genes), (genes, cells)), shape=(self.n_genes, self.n_cells), dtype=np.uint16)
		self.total_umis = np.array(self.matrix.sum(axis=0))[0]
		return self.matrix

	def remove_empty_beads(self, expected_n_cells: int) -> None:
		logging.info("Calling cells")
		self.ambient_umis, self.ambient_pvalue = call_cells(self.matrix.tocsc(), expected_n_cells)
		self.valid_cells = (self.ambient_pvalue < 0.01) | (self.total_umis > 1500)
		self.matrix = self.matrix.tocsr()[:, self.valid_cells]
		self.cell_ids = self.cell_ids[self.valid_cells]
		self.n_cells = self.valid_cells.sum()
		for name in self.layers.keys():
			self.layers[name] = self.layers[name].tocsc()[:, self.valid_cells]
		logging.info(f"Found {self.n_cells} valid cells and ~{int(self.ambient_umis)} ambient UMIs.")

	def count_layer(self, layer_name: str, layer_fragments_file: str) -> sparse.coo_matrix:
		fragments_idx = {f: i for i, f in enumerate(self.fragments)}
		with open(layer_fragments_file) as f:
			include_fragments = set([fragments_idx[x[:-1]] for x in f.readlines()])
		logging.info(f"Counting pseudoalignments for layer '{layer_name}'")
		# Figure out which of the equivalence classes are relevant
		include_ec = {}
		for ec, tids in self.equivalence_classes.items():
			if any([tid in include_fragments for tid in tids]):
				include_ec[ec] = True
			else:
				include_ec[ec] = False
		m = IncrementalSparseMatrixUInt16((self.n_genes, self.n_cells))
		for ix in range(self.n_records):
			if not self.bus_valid[ix]:
				continue
			gene = self.bus_gene[ix]
			if include_ec[self.bus["equivalence_class"][ix]]:
				m.append(gene, self.cell_for_barcode_idx[self.bus["barcode"][ix]], 1)
		self.layers[layer_name] = m.tocoo()
		return self.layers[layer_name]

	def save(self, out_file: str, sample_id: str, samples_metadata_file: str) -> None:
		logging.info("Saving")
		row_attrs = {}
		# Transpose the gene metadata
		for i, attr in enumerate(self.gene_metadata_attributes):
			row_attrs[attr] = np.array([v[i] for v in self.genes.values()])

		# Create cell attributes
		col_attrs = {
			"CellID": np.array([sample_id + "_" + twobit_to_dna(int(cid), 16) for cid in self.cell_ids]),
			"TotalUMIs": self.total_umis[self.valid_cells]
		}

		# Load sample metadata
		metadata = load_sample_metadata(samples_metadata_file, sample_id)
		global_attrs = {
			**{
				"SampleID": sample_id,
				"AmbientUMIs": self.ambient_umis,
				"RedundantReadFraction": 1 - self.bus_valid.sum() / self.n_records,
				"AmbientPValue": self.ambient_pvalue,
				"BarcodeTotalUMIs": self.total_umis,
				"CellBarcodes": self.valid_cells
			}, **metadata}

		layers = self.layers.copy()
		layers[""] = self.matrix
		create(out_file, layers, row_attrs, col_attrs, file_attrs=global_attrs)


def execute(cmd: List[str], synchronous: bool = False) -> Generator:
	if synchronous:
		yield os.popen(" ".join(cmd)).read()
	else:
		popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)  # type: ignore
		for stdout_line in iter(popen.stdout.readline, ""):
			yield stdout_line
		popen.stdout.close()
		return_code = popen.wait()
		if return_code:
			raise subprocess.CalledProcessError(return_code, cmd)


def create_from_fastq(out_file: str, sample_id: str, fastqs: List[str], index_path: str, samples_metadata_file: str, n_threads: int = 1, temp_folder: str = None, synchronous: bool = False) -> None:
	"""
	Args:
		technology			String like "10xv2" or None to read the technology from the sample metadata file
		expected_n_cells	Expected number of cells captured in the sample, or None to read the number from the sample metadata file
		samples_metadata_file	Path to tab-delimited file with one header row OR path to sqlite database with one table called "sample"
	
	Remarks:
		Samples metadata table should contain these columns:
			name			Sample name (i.e. sample id)
			chemistry		10x chemistry version (v1, v2 or v3)
			targetnumcells	Number of cells expected in the sample
	"""
	manifest_file = os.path.join(index_path, "manifest.json")
	if not os.path.exists(manifest_file):
		raise ValueError(f"Manifest file 'manifest.json' was missing from index at '{index_path}'")
	for fastq in fastqs:
		if not os.path.exists(fastq):
			raise ValueError(f"Fastq file '{fastq}' was not found")
	if not os.path.exists(samples_metadata_file):
		raise ValueError("Samples metadata file not found")
	with open(manifest_file) as f:
		manifest = json.load(f)
	
	metadata = load_sample_metadata(samples_metadata_file, sample_id)

	if (("technology" not in metadata) and ("chemistry" not in metadata)) or ("targetnumcells" not in metadata):
		print(metadata.keys())
		raise ValueError("Samples metadata must contain columns 'targetnumcells' and either 'chemistry' or 'technology'")
	if "technology" in metadata:
		technology = metadata["technology"]
	else:
		technology = "10x" + metadata["chemistry"]
	try:
		expected_n_cells = int(metadata["targetnumcells"])
	except:
		expected_n_cells = 5000

	whitelist_file: Optional[str] = os.path.join(index_path, f"{technology}_whitelist.txt")
	if not os.path.exists(whitelist_file):  # type: ignore
		logging.warning(f"Barcode whitelist file {whitelist_file} not found in index folder at '{index_path}'; barcode correction will be skipped.")
		whitelist_file = None

	with TemporaryDirectory() as d:
		if temp_folder is not None:
			d = temp_folder
			if not os.path.exists(d):
				os.mkdir(d)
		cmd = ["kallisto", "bus", "-i", os.path.join(index_path, manifest["index_file"]), "-o", d, "-x", technology, "-t", str(n_threads)] + fastqs
		logging.info(" ".join(cmd))
		for line in execute(cmd, synchronous):
			if line != "\n":
				logging.info(line[:-1])

		run_info: Optional[Dict[str, str]] = None
		try:
			with open(os.path.join(d, "run_info.json")) as f:
				run_info = json.load(f)
		except json.JSONDecodeError as e:
			with open(os.path.join(d, "run_info.json")) as f:
				for line in f:
					logging.error(line)
			logging.error(f"Error decoding run_info.json: {e}")
		bus = BusFile(
			os.path.join(d, "output.bus"),
			os.path.join(index_path, manifest["gene_metadata_file"]),
			manifest["gene_metadata_key"],
			os.path.join(index_path, manifest["fragments_to_genes_file"]),
			os.path.join(d, "matrix.ec"),
			os.path.join(d, "transcripts.txt")
		)
		logging.info(f"Found {bus.n_records:,} records for {bus.n_genes:,} genes and {bus.n_cells:,} uncorrected cell barcodes.")
		if whitelist_file is None:
			logging.warning("Not correcting barcodes, because whitelist file was not provided.")
		else:
			logging.info("Correcting cell barcodes")
		bus.correct(whitelist_file)
		logging.info(f"Found {bus.n_cells:,} corrected cell barcodes.")
		logging.info("Removing redundant reads using UMIs")
		bus.deduplicate()
		seq_sat = 1 - bus.bus_valid.sum() / bus.n_records
		logging.info(f"{int(seq_sat * 100)}% sequencing saturation.")
		bus.count()
		logging.info(f"Found {bus.matrix.count_nonzero():,} UMIs.")
		for layer, layer_def in manifest["layers"].items():
			bus.count_layer(layer, os.path.join(index_path, layer_def))
			logging.info(f"Found {bus.layers[layer].count_nonzero():,} UMIs.")
		bus.remove_empty_beads(expected_n_cells)
		logging.info(f"Creating loom file '{out_file}'")
		bus.save(out_file, sample_id, samples_metadata_file)
		with connect(out_file) as ds:
			for ra in ds.row_attrs:  # For some reason it happens that e.g. ds.ra.Alias == None, which fails downstream work.
				if ds.rowattrs[ra] is None:
					del ds.rowattrs[ra]
			ds.attrs.Species = manifest["species"]
			ds.attrs.Saturation = seq_sat
			if run_info is not None:
				ds.attrs.NumReadsProcessed = int(run_info["n_processed"])
				ds.attrs.NumPseudoaligned = int(run_info["n_pseudoaligned"])
				ds.attrs.KallistoCommand = run_info["call"]
				ds.attrs.KallistoVersion = run_info["kallisto_version"]

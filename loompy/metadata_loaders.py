from typing import Dict, Union, Iterable
import os.path, re
import sqlite3 as sqlite
import numpy as np

def load_gene_metadata(gtf_file : str) -> Dict[str, Dict[str, Union[int, str]]]:
	"""
	Read gene metadata from a GTF file.

	Args:
	  gtf_file (str):             path to GTF file

	Returns:
	  A Dict with each Accession (gtf "gene_id") pointing to a Dict of metadata keys -> values
	"""
	if not os.path.exists(gtf_file):
		raise ValueError(f"Gene metadata file '{gtf_file}' not found.")
	regex_genetype = re.compile('gene_biotype "([^"]+)"')
	regex_geneid = re.compile('gene_id "([^"]+)"')
	regex_genename = re.compile('gene_name "([^"]+)"')
	geneid2annots = {}
	for line in open(gtf_file).readlines():
		if line.startswith('#'):
			continue
		fields = line.rstrip().split('\t')
		chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
		if feature_type == "gene":
			genename = geneid = regex_geneid.search(tags).group(1)
			_genename_search = regex_genename.search(tags)
			if _genename_search:
				genename = _genename_search.group(1)
			_genetype_search = regex_genetype.search(tags)
			genetype = _genetype_search.group(1) if _genetype_search else "n/a"
			chrid, start, end = fields[0], int(fields[3]), int(fields[4])
			attrs = { "Gene": genename, "Accession": geneid, "Biotype": genetype, \
                                                  "Chromosome": chrid, "Start": start, "End": end }
			geneid2annots[geneid] = attrs
			geneid2annots[genename] = attrs
	return geneid2annots

def make_row_attrs_from_gene_annotations(acc2annots : Dict[str, Dict[str, Union[int, str]]], ordered_features : Iterable[str]) -> Dict[str, np.ndarray]:
	"""
         Construct loom row attributes corresponding to ordered_features from load_gene_metadata output.

	Args:
	  acc2annots (Dict):          output from load_gene_metadata
          ordered_features (str):     the accessions (should match the gtf "gene_id") in matrix row order

	Returns:
          A row attribute dictionary of attr->numpy arrays to assign to a Loom object.
	"""
	ra = {}
	first_annot = next(iter(acc2annots.values()))
	ra_attrs = list(first_annot.keys())
	n_genes = len(ordered_features)
	for ra_attr in ra_attrs:
		ra[ra_attr] = np.zeros((n_genes,), dtype = object)
	for idx, geneid in enumerate(ordered_features):
		try:
			annots = acc2annots[geneid]
		except KeyError:
			annots = { "Gene": geneid, "Accession": geneid, "Biotype": "n/a", "Chromosome": "Un", "Start": "0", "End": "0" }
		for ra_attr in ra_attrs:
			ra[ra_attr][idx] = annots[ra_attr]
	return ra

def make_row_attrs_from_gene_metadata(gtf_file : str, ordered_features : Iterable[str]) -> Dict[str, np.ndarray]:
	"""
        Read gene metadata from a GTF file and construct loom row attributes corresponding to ordered_features.

	Args:
	  gtf_file (str):             path to GTF file
          ordered_features (str):     the accessions (should match the gtf "gene_id") in matrix row order

	Returns:
          A row attribute object ready to assign to a Loom object.
	"""
	acc2annots = load_gene_metadata(gtf_file)
	return make_row_attrs_from_gene_annotations(acc2annots, ordered_features)

def load_sample_metadata(path: str, sample_id: str) -> Dict[str, str]:
        """
        Read sample metadata from either an sqlite '.db' database or a tab-delimited file.
        The tab-file has to have  one sample per line and a header with a column
        labelled 'Name' or 'SampleId' in which a match to function argument sample_id should be found.

        Args:
          path (str):             path to sqlite database (.db) or tab-delimited metadata file
          sample_id (str):        Id of sample to annotate

        Returns:
              A Dict of metadata keys -> values
        """
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


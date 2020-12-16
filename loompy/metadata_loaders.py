from typing import Dict, Union
import os.path, re
import sqlite3 as sqlite

def load_gene_metadata(gtf_file : str) -> Dict[str, Dict[str, Union[int, str]]]:
        """
        Read gene metadata from a GTF file.
        Args:
          indir (str):             path to GTF file
        Returns:
              A Dict with each GeneId pointing to a Dict of metadata keys -> values
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
        		genetype = regex_genetype.search(tags).group(1)
        		chrid, start, end = fields[0], int(fields[3]), int(fields[4])
        		geneid2annots[geneid] = { "Gene:": genename, "Accession": geneid, "Biotype": genetype, \
        		                          "Chromosome": chrid, "Start": start, "End": end }
        return geneid2annots

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


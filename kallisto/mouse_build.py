from collections import defaultdict

d = "./"

mgiID2MRK_ENSE = {}
enseID2mgiID = {}
with open(d + "inputs/MRK_ENSEMBL.rpt") as f:
	for line in f:
		items = line[:-1].split("\t")
		mgiID = items[0]
		mgiID2MRK_ENSE[mgiID] = items
		enseID = items[5]
		enseID2mgiID[enseID] = mgiID

enseID2MGI_GMC = {}
with open(d + "inputs/MGI_Gene_Model_Coord.rpt") as f:
	MGI_GMC_headers = f.readline()[:-1].split("\t")
	for line in f:
		items = line[:-1].split("\t")
		enseID = items[10]
		enseID2MGI_GMC[enseID] = items

mgiID2MRK_Seq = {}
with open(d + "inputs/MRK_Sequence.rpt") as f:
        MRK_Seq_headers = f.readline()[:-1].split("\t")
        for line in f:
                items = line[:-1].split("\t")
                mgiID = items[0]
                mgiID2MRK_Seq[mgiID] = items

enseID2mart = {}
with open(d + "inputs/mart_export.txt") as f:
	mart_headers = f.readline()[:-1].split("\t")
	for line in f:
		items = line[:-1].split("\t")
		enseID = items[0]
		enseID2mart[enseID] = items

geneSymbol2TF = {}
with open(d + "inputs/TF_TcoF-DB.tsv") as f:
	TF_headers = f.readline()[:-1].split("\t")
	for line in f:
		items = line[:-1].split("\t")
		geneSymbol = items[0]
		geneSymbol2TF[geneSymbol] = items

geneSymbol2Regulated = defaultdict(list)
with open(d + "inputs/trrust_rawdata.mouse.tsv") as f:
	for line in f:
		items = line[:-1].split("\t")
		TFSymbol = items[0]
		geneSymbol2Regulated[TFSymbol].append(items[1])

with open(d + "gencode.vM23.metadata.tab", "w") as fout:
	fout.write("\t".join([
		"Accession",
		"AccessionVersion",
		"Gene",
		"FullName",
		"GeneType",
		"HgncID",
		"Chromosome",
		"Strand",
		"ChromosomeStart",
		"ChromosomeEnd",
		"LocusGroup",
		"LocusType",
		"Location",
		"LocationSortable",
		"Aliases",
		"VegaID",
		"UcscID",
		"RefseqID",
		"CcdsID",
		"UniprotID",
		"PubmedID",
		"MgdID",
		"RgdID",
		"CosmicID",
		"OmimID",
		"MirBaseID",
		"IsTFi (TcoF-DB)",
		"DnaBindingDomain",
		"Regulates (TRRUST)"
	]))
	fout.write("\n")
	with open(d + "inputs/gencode.vM23.primary_assembly.annotation.gtf") as f:
		for line in f:
			if line.startswith("##"):
				continue
			items = line[:-1].split("\t")
			if items[2] != "gene":
				continue
			extra = {x.strip().split(" ")[0]: x.strip().split(" ")[1].strip('"') for x in items[8].split(";")[:-1]}
			enseID = extra["gene_id"].split(".")[0]
			geneSymbol = extra.get("gene_name", "")
			fout.write("\t".join([
				enseID,
				extra["gene_id"],
				geneSymbol,
				enseID2MGI_GMC[enseID][3] if enseID in enseID2MGI_GMC else "",  # full name
				extra["gene_type"],  # gene type from gencode
				"",  # HGNC id
				items[0],  # Chromosome
				items[6],
				items[3],  # Start
				items[4],  # End
				"",  # Locus group
				mgiID2MRK_ENSE[mgiID][8],  # Locus type
				"",  # Location
				"",  # Location, sortable
				"",  # Aliases
				enseID2mart[enseID][5] if enseID in enseID2mart else "",  # VEGA id
				enseID2mart[enseID][4] if enseID in enseID2mart else "",  # UCSC id
				mgiID2MRK_Seq[mgiID][12],  # Refseq id
				enseID2mart[enseID][6] if enseID in enseID2mart else "",  # CCDS id
				mgiID2MRK_Seq[mgiID][14],  # Uniprot id
				"",  # Pubmed id
				"",  # MGD id
				"",  # RGD id
				"",  # COSMIC id
				"",  # OMIM id
				"",  # MIRbase id
				"True" if (geneSymbol in geneSymbol2TF) else "False", #  IsTF?
				"", # DBD
				",".join(geneSymbol2Regulated[geneSymbol])  #  TF regulated genes
				]))
			fout.write("\n")

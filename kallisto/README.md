# Building an annotated kallisto index

For the human genome see the notebook subdirectory.

## Building the mouse kallisto index

These instructions work on Linux (tested on CentOS7).

1. Make sure packages bedtools and kallisto are installed on the system.
  - bedtools from https://bedtools.readthedocs.io/en/latest/content/installation.html
  - kallisto from https://pachterlab.github.io/kallisto/download.html

2. Create your working directory, 'cd' there, and put the files above there.

3. Download and preprocess input files:

`bash mouse_download.sh`

  This will create a directory "inputs" and put some files there as well as in the current directory.

4. Download "BrowseTF  TcoF-DB.xlsx" from https://tools.sschmeier.com/tcof/browse/?type=tcof&species=mouse&class=all# by clicking the "Excel" button. (Main page is https://tools.sschmeier.com/tcof/home/).
  Open the file in Excel and save tab-separated as "inputs/TcoF-DB.tsv".

5. You need to download some annotations for Mouse GRCm38 from BioMart (https://m.ensembl.org/biomart) Open this link in a new browser tab:

http://www.ensembl.org/biomart/martview/7c9b283e3eca26cb81449ec518f4fc14?VIRTUALSCHEMANAME=default&ATTRIBUTES=mmusculus_gene_ensembl.default.feature_page.ensembl_gene_id|mmusculus_gene_ensembl.default.feature_page.ensembl_gene_id_version|mmusculus_gene_ensembl.default.feature_page.ensembl_transcript_id|mmusculus_gene_ensembl.default.feature_page.ensembl_transcript_id_version|mmusculus_gene_ensembl.default.feature_page.ucsc|mmusculus_gene_ensembl.default.feature_page.vega_translation|mmusculus_gene_ensembl.default.feature_page.ccds&FILTERS=&VISIBLEPANEL=resultspanel

  On this BioMart page, click the "Go" button, and save the downloaded "mart_export.txt" file as "inputs/mart_export.txt".
  The file should contain the following columns in the header:
      Gene stable ID        Gene stable ID version  Transcript stable ID    Transcript stable ID version    UCSC Stable ID  Vega translation ID     CCDS ID

  If the link fails, you need to manually select the proper dataset and columns from the https://m.ensembl.org/biomart webpage and download:
  *  Select Dataset "Ensembl Genes 101"/"Mouse genes GRCm38".
  *  Select Attributes as in columns above: First 4 should be auto-selected. Select the following 3 from the "EXTERNAL" section, clicking the 3 boxes in the order above.
  * Click "Results", export using "Go" button, and save to "inputs/mart_export.txt". 

6. Run the annotation assembly script:

`python mouse_build.py`

7. Create the "manifest.json" file or use the one supplied above. It should contain:
```
{
    "species": "Mus musculus",
    "index_file": "gencode.vM23.fragments.idx",
    "gene_metadata_file": "gencode.vM23.metadata.tab",
    "gene_metadata_key": "AccessionVersion",
    "fragments_to_genes_file": "fragments2genes.txt",
    "layers": {
        "unspliced": "unspliced_fragments.txt",
        "spliced": "spliced_fragments.txt"
    }
}
```

8. Run the fragment generator script:

`python mouse_generate_fragments.py`

9. Build the kallisto index:

`kallisto index -i gencode.vM23.fragments.idx -k 31 inputs/gencode.vM23.fragments.fa`

10. Refer to the notebook for human for more info on the output.

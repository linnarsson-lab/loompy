mkdir inputs

# Make sure these manual steps have been done:
# Download "BrowseTF  TcoF-DB.xlsx" from https://tools.sschmeier.com/tcof/browse/?type=tcof&species=mouse&class=all# (a button at https://tools.sschmeier.com/tcof/home/)
#   Open the file in Excel and save tab-separated as "inputs/TcoF-Db.tsv"
#
# You need to import data from BioMart using this link:
# http://www.ensembl.org/biomart/martview/7c9b283e3eca26cb81449ec518f4fc14?VIRTUALSCHEMANAME=default&ATTRIBUTES=mmusculus_gene_ensembl.default.feature_page.ensembl_gene_id|mmusculus_gene_ensembl.default.feature_page.ensembl_gene_id_version|mmusculus_gene_ensembl.default.feature_page.ensembl_transcript_id|mmusculus_gene_ensembl.default.feature_page.ensembl_transcript_id_version|mmusculus_gene_ensembl.default.feature_page.ucsc|mmusculus_gene_ensembl.default.feature_page.vega_translation|mmusculus_gene_ensembl.default.feature_page.ccds&FILTERS=&VISIBLEPANEL=resultspanel
# by clicking "Go" button, and saving the downloaded "mart_export.txt" file in "inputs/mart_export.txt".
# The file should contain the following columns:
# Gene stable ID	Gene stable ID version	Transcript stable ID	Transcript stable ID version	UCSC Stable ID	Vega translation ID	CCDS ID

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.primary_assembly.annotation.gtf.gz
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-april-2014_rc.txt
wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt

zcat gencode.vM23.primary_assembly.annotation.gtf.gz | gawk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,0,$7}}' | tr -d '";' > gencode.vM23.primary_assembly.annotation.bed
bedtools sort -i gencode.vM23.primary_assembly.annotation.bed > gencode.vM23.primary_assembly.annotation.sorted.bed
bedtools merge -i gencode.vM23.primary_assembly.annotation.sorted.bed -s -c 4 -o collapse > gencode.vM23.primary_assembly.annotation.merged.bed
gunzip GRCm38.primary_assembly.genome.fa.gz 
bedtools getfasta -name -fo gencode.vM23.unspliced.fa -fi GRCm38.primary_assembly.genome.fa -bed gencode.vM23.primary_assembly.annotation.sorted.bed

mv 737K-april-2014_rc.txt 10xv1_whitelist.txt
mv 737K-august-2016.txt 10xv2_whitelist.txt
gunzip 3M-february-2018.txt.gz 
mv 3M-february-2018.txt 10xv3_whitelist.txt 

mv GRCm38.primary_assembly.genome.fa* inputs/
mv gencode.vM23.unspliced.fa inputs/
mv gencode.vM23.primary_assembly.annotation.bed inputs/
mv gencode.vM23.primary_assembly.annotation.sorted.bed inputs/
gunzip gencode.vM23.primary_assembly.annotation.gtf.gz 
mv gencode.vM23.primary_assembly.annotation.gtf inputs/

cd inputs
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
wget https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv
gunzip gencode.vM23.transcripts.fa.gz 
wget http://www.informatics.jax.org/downloads/reports/MGI_Gene_Model_Coord.rpt
wget http://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt
wget http://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt
cd ..


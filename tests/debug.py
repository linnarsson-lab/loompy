import sys
import logging
import loompy

# Create logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Create STDERR handler
handler = logging.StreamHandler(sys.stderr)
# ch.setLevel(logging.DEBUG)

# Create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)

# Set STDERR handler as the only handler 
logger.handlers = [handler]

d = "/Users/stelin/kallisto_GRCh38/"

fastqs = [
	d + "HLFGJBCXY/10X_17_029_S2_L002_R1_001.fastq.gz",
	d + "HLFGJBCXY/10X_17_029_S2_L002_R2_001.fastq.gz"
#	d + "HL73JBCXY/10X_17_029_S2_L002_R1_001.fastq.gz",
#	d + "HL73JBCXY/10X_17_029_S2_L002_R2_001.fastq.gz"
]
#loompy.create_from_fastq(d + "10X_17_029.loom", "10X_17_029", fastqs, d + "human_GRCh38_gencode.v31", d + "samples.tab", n_threads=6)
loompy.create_from_fastq(d + "10X_17_029.loom", "10X_17_029", fastqs, d + "human_GRCh38_gencode.v31", "/Users/stelin/sqlite3_chromium.db", n_threads=6)

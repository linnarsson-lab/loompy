import sys, os

extent = 600  # how many bases away from polya to include
min_len = 90  # how many non-repeat bases required to make a transcript

from typing import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def find_polys(seq: SeqRecord, c: str = "A", n: int = 15) -> List[Tuple[int, int]]:
	found = []
	count = seq[:n].count(c)  # Count occurences in the first k-mer
	if count >= n - 1:  # We have a match
		found.append(0)
	ix = 0
	while ix < len(seq) - n - 1:
		if seq[ix] == c:  # Outgoing base
			count -= 1
		if seq[ix + n] == c:  # Incoming base
			count += 1
		ix += 1
		if count >= n - 1:  # We have a match
			found.append(ix)
	
	sorted_by_lower_bound = [(f, f + n) for f in found]
	# merge intervals (https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals)
	merged = []
	for higher in sorted_by_lower_bound:
		if not merged:
			merged.append(higher)
		else:
			lower = merged[-1]
			# test for intersection between lower and higher:
			# we know via sorting that lower[0] <= higher[0]
			if higher[0] <= lower[1]:
				upper_bound = max(lower[1], higher[1])
				merged[-1] = (lower[0], upper_bound)  # replace by merged interval
			else:
				merged.append(higher)
	return merged

polyAs = {}
polyTs = {}
for fasta in SeqIO.parse(open("inputs/gencode.vM23.unspliced.fa"),'fasta'):
	gene_id = fasta.id
	intervals = find_polys(fasta.seq, c="A", n=14)
	if len(intervals) > 0:
		polyAs[gene_id] = intervals
	# Collect fragments on the opposite strand, downstream of poly-Ts (not sure if such reads really happen?)
	intervals = find_polys(fasta.seq, c="T", n=14)
	if len(intervals) > 0:
		polyTs[gene_id] = intervals

tr2g = {}
with open("inputs/gencode.vM23.primary_assembly.annotation.gtf") as f:
	for line in f:
		if "\ttranscript\t" in line:
			items = line.split("; ")
			chrom, _, _, start, end, _, strand, _, gid = items[0].split("\t")
			gene_id = gid.split('"')[1]
			transcript_id = items[1].split('"')[1]
			gene_type = items[2].split('"')[1]
			gene_name = items[3].split('"')[1]
			tr2g[transcript_id] = (chrom, start, end, strand, gene_id, gene_type, gene_name)

count = 0
with open("fragments2genes.txt", "w") as ftr2g:
	with open("inputs/gencode.vM23.fragments.fa", "w") as fout:
		# Write the nascent fragments, with one partial transcript per internal poly-A/T site
		with open("unspliced_fragments.txt", "w") as fucapture:
			for fasta in SeqIO.parse(open("inputs/gencode.vM23.unspliced.fa"),'fasta'):  # Note we're in the masked file now
				gene_id = fasta.id
				if gene_id in polyAs:
					for interval in polyAs[gene_id]:
						seq = str(fasta.seq[max(0, interval[0] - extent):interval[0]])
						#seq = seq.translate(tr).strip("N")
						if len(seq) >= min_len:
							count += 1
							transcript_id = f"{gene_id}.A{interval[0]}"
							trseq = SeqRecord(Seq(seq), transcript_id, '', '')
							fout.write(trseq.format("fasta"))
							ftr2g.write(f"{transcript_id}\t{gene_id}\n")
							fucapture.write(f"{transcript_id}\n")
				if gene_id in polyTs:
					for interval in polyTs[gene_id]:
						seq = str(fasta.seq[interval[1]:interval[1] + extent].reverse_complement())
						#seq = seq.translate(tr).strip("N")
						if len(seq) >= min_len:
							count += 1
							transcript_id = f"{gene_id}.T{interval[0]}"
							trseq = SeqRecord(Seq(seq), transcript_id, '', '')
							fout.write(trseq.format("fasta"))
							ftr2g.write(f"{transcript_id}\t{gene_id}\n")
							fucapture.write(f"{transcript_id}\n")
		# Write the mature fragments, covering the 3' end of each mature transcript
		with open("spliced_fragments.txt", "w") as fscapture:
			for fasta in SeqIO.parse(open("inputs/gencode.vM23.transcripts.fa"),'fasta'):  # Note we're in the masked file now
				transcript_id = fasta.id.split("|")[0]
				gene_id = fasta.id.split("|")[1]
				attrs = tr2g[transcript_id]
				seq = str(fasta.seq[-extent:])
				if len(seq) >= min_len:
					count += 1
					trseq = SeqRecord(Seq(seq), f"{transcript_id}.{count} gene_id:{attrs[4]} gene_name:{attrs[6]}", '', '')
					fout.write(trseq.format("fasta"))
					ftr2g.write(f"{transcript_id}.{count}\t{attrs[4]}\n")
					fscapture.write(f"{transcript_id}.{count}\n")


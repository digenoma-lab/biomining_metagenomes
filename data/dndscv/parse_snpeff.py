#!/usr/bin/env python3
"""
Extract mutations from the VCF file obtained with SnpEff,
generating a parsed file for dNdScv package as output

python3 parse_snpeff.py S15.ann.vcf > S15_mutations.txt
"""

import sys

# Read VCF file from path
vcf_name = sys.argv[1]
vcf_file = open(vcf_name, "r")

i = 0
chrom = []
pos = []
ref = []
mut = []
sid = []

# Parse annotations
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S15_sorted.bam
for line in vcf_file:
	line = line.strip()
	if line[:2] != "##" and line[:6] != "#CHROM":
		#i += 1
		data = line.split("\t")
		#keep attributes
		chrom.append(data[0])
		pos.append(data[1])
		ref.append(data[3])
		mut.append(data[4])
		sid.append("S15") #data[2]
		#info = data[7].split(";")
		#print(len(info))
		#last ; is ''ANN='', then split by '|'
		#example: ANN=C|synonymous_variant|LOW|
		#annot = info[-1].split("|")
		#if info[0][:4] == "ANN=":

vcf_file.close()

# Print header and output columns
print("sampleID\tchr\tpos\tref\tmut")
for x in range(len(chrom)):
	print(sid[x] + "\t" + chrom[x] + "\t" + pos[x] + "\t" + ref[x] + "\t" + mut[x])

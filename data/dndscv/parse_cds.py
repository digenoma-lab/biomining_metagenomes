#!/usr/bin/env python3
"""
Take as input a GFF3 file (ie. from Bakta) and parse the annotation of genes,
printing a tab-delimited CDS table with defined columns as output.

python3 parse_cds.py data/S15_bins/genes.gff > cds_table.txt
"""

import sys

# Read GFF file from path
gff_name = sys.argv[1]
gff_file = open(gff_name, "r")

# Lists for output..
i = 0
j = 0
chrom = []
chrom_start = []
chrom_stop = []
cds_start = []
#cds_stop = []
length = []
strand = []
gene_id = []
gene_name = []

# Read CDS annotation: genes, chrom (bins), positions, etc.
# Since dNdScv is focused on eukaryotic genomes, we need to make adjustments for bacteria (no splicing)
#https://github.com/im3sanger/dndscv/issues/62
for line in gff_file:
	line = line.strip()
	#NODE_1370_length_37472_cov_9.854745	Prodigal:002006	CDS	1078	3408	.	+	0	ID=MBCHOMPH_01850;inference=ab initio prediction:Prodigal:002006;locus_tag=MBCHOMPH_01850;product=hypothetical protein
	if line[:5] == "NODE_":
		#i += 1
		data = line.split("\t")
		# Only features with type=CDS
		if data[2] == "CDS":
			#j += 1
			chrom.append(data[0]) #chrom name
			chrom_start.append(data[3]) #chrom start
			chrom_stop.append(data[4]) #chrom stop
			cds_start.append("1") #CDS start
			#pos = data[0].split("_") #get max pos from range
			#cds_stop.append(pos[3]) #chrom stop (from gff ##sequence-region, last col)
			l = int(data[4]) - int(data[3]) + 1 #cast str to int
			length.append(l)
			if data[6] == "+": #strand + or -
				strand.append("1")
			else:
				strand.append("-1")
			att = data[8].split(";") #parse to ID=gene.name; Name=name;product=hypothetical|description
			gene_id.append(att[0][3:])
			ini = [i for i in range(len(att)) if att[i].startswith("Name=")] #get index
			if ini:
				gene_name.append(att[ini[0]][5:])
			else:
				gene_name.append(att[0][3:])

#print(len(chrom))
gff_file.close()

# Print header with required columns
print("gene.id\tgene.name\tcds.id\tchr\tchr.coding.start\tchr.coding.end\tcds.start\tcds.end\tlength\tstrand")
for x in range(len(chrom)):
	print(gene_id[x] + "\t" + gene_id[x] + "\t" + gene_id[x] + "\t" + chrom[x] + "\t" + chrom_start[x] + "\t" + chrom_stop[x] + "\t" + cds_start[x] + "\t" + str(length[x]) + "\t" + str(length[x]) + "\t" + strand[x])

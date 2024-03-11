#!/usr/bin/env Rscript

#---------------------------------------------------------
# Systemix MAGs: Get the QC from a list of -refined- bins
# 
# Moises A. Rojas 06/11/2023 - mod. 06/03/2024
#---------------------------------------------------------

# Set working dir
setwd("/mnt/beegfs/home/mrojas/DiGenomaLab/Systemix/analysis/mags-all/analysis/2-statistics/")

# Read bin list
bins <- read.csv("bins_refined_list.tsv", sep = "\t", stringsAsFactors = FALSE)
selected_bins <- bins$bin

# QC parameters to keep
completeness  <- vector()
contamination <- vector()
strain_het    <- vector()
genome_size   <- vector()
contigs       <- vector()
n50_contigs   <- vector()
gc            <- vector()
pred_genes    <- vector()

# Read the QC file from CheckM
checkm_qc <- read.csv("checkm_qa.tsv", sep = "\t", stringsAsFactors = FALSE)

# Save attributes for each bin
checkm_bin <- checkm_qc$bin
checkm_com <- checkm_qc$completeness
checkm_con <- checkm_qc$contamination
checkm_st  <- checkm_qc$strain_heterogeneity
checkm_gs  <- checkm_qc$genome_size
checkm_nco <- checkm_qc$contigs
checkm_n50 <- checkm_qc$N50_contigs
checkm_gc  <- checkm_qc$gc
checkm_pg  <- checkm_qc$predicted_genes

# Match positions for each bin, saving attributes
for (i in 1:length(selected_bins)) { #428
	for (j in 1:length(checkm_bin)) { #1840
		if (selected_bins[i] == checkm_bin[j]) {
			#print(j) #pos
			completeness  <- append(completeness,  checkm_com[j])
			contamination <- append(contamination, checkm_con[j])
			strain_het    <- append(strain_het,    checkm_st[j])
			genome_size   <- append(genome_size,   checkm_gs[j])
			contigs       <- append(contigs,       checkm_nco[j])
			n50_contigs   <- append(n50_contigs,   checkm_n50[j])
			gc            <- append(gc,            checkm_gc[j])
			pred_genes    <- append(pred_genes,    checkm_pg[j])
		}
	}
}

# Create output dataframe and save to file
df <- data.frame(bin = selected_bins,
				completeness = completeness,
				contamination = contamination,
				strain_heterogeneity = strain_het,
				genome_size = genome_size,
				contigs = contigs,
				N50 = n50_contigs,
				gc = gc,
				predicted_genes = pred_genes)
write.table(df, "checkm_parsed_bins.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

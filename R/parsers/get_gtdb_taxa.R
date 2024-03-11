#!/usr/bin/env Rscript

#-----------------------------------------------
# Populate the GTDB-Tk of MAGs from a reference
# Moises A. Rojas 13/10/2023
#-----------------------------------------------

setwd("/mnt/beegfs/home/mrojas/DiGenomaLab/Systemix/analysis/mags-all/analysis/1-taxonomy")

# Read the list of bins (last nf-core/mag run: mags-all/)
bins_all <- read.csv("bins_refined_list.tsv", sep = "\t", stringsAsFactors = FALSE)
head(bins_all)

bins_genome <- bins_all$bin

# Read the GTDB-Tk classification file
gtdb_file <- read.csv("gtdbtk_summary.tsv", sep = "\t", stringsAsFactors = FALSE)
head(gtdb_file)

# Save attributes
gtdb_genome <- gtdb_file$user_genome
gtdb_taxa   <- gtdb_file$classification

# Set output attributes
bin     <- vector()
taxa    <- vector()
kingdom <- vector()
phylum  <- vector()
class   <- vector()
order   <- vector()
family  <- vector()
genus   <- vector()
species <- vector()

# Parse taxonomy array:
# d__Bacteria;p__Nitrospirota_A;c__Leptospirillia;o__Leptospirillales;f__Leptospirillaceae;g__;s__
taxa_split <- strsplit(gtdb_taxa, split = ";", fixed = TRUE)

# Match positions for each bin
# Then, save taxa attributes
for (i in 1:length(bins_genome)) { # 428
	for (j in 1:length(gtdb_genome)) { # 377
		if (bins_genome[i] == gtdb_genome[j]) {
			#print(j) #pos
			bin     <- append(bin,     gtdb_genome[[j]])
			taxa    <- append(taxa,    gtdb_taxa[[j]])
			kingdom <- append(kingdom, taxa_split[[j]][1])
			phylum  <- append(phylum,  taxa_split[[j]][2])
			class   <- append(class,   taxa_split[[j]][3])
			order   <- append(order,   taxa_split[[j]][4])
			family  <- append(family,  taxa_split[[j]][5])
			genus   <- append(genus,   taxa_split[[j]][6])
			species <- append(species, taxa_split[[j]][7])
		}
	}
}

# Create output dataframe and save to file
df <- data.frame(genome  = bin,
				taxonomy = taxa,
				kingdom  = kingdom,
				phylum   = phylum,
				class    = class,
				order    = order,
				family   = family,
				genus    = genus,
				species  = species)

write.table(df, "gtdbtk_classification.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

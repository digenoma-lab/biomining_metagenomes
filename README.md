# biomining_metagenomes

Repository that holds data, scripts and figures regarding to the mining MAGs (Cauquenes tailing) article: **"Genome-resolved metagenomics and evolutionary analysis reveal conserved metabolic adaptations in extremophile communities from a copper mining tailing"**.

## Section 1: MAG taxonomy, QC and abundances

We classified a set of 44 bacterial metagenome-assembled genomes (MAGs) obtained from the ```S15``` sample (intermediate tailing) using [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk). We also assessed the quality metrics (generated with [CheckM](https://github.com/Ecogenomics/CheckM)), and the abundances/coverage of the mining bins.

See [main figure #1](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_1.png).

Also [supplementary figure #1](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s1.png),
[supplementary figure #2](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s2.png).

## Section 2: Functional annotation

Next, we compared the 44 mining MAGs with a subset of MAGs from nine different non-mining ecosystems ([SMAG catalog](https://microbma.github.io/project/SMAG.html), bacterial bins) to assess their capacities for selected mining pathways (copper, iron, and sulfur proteins).

See [main figure #2](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_2.png).

Also [supplementary figure #3](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s3.png),
[supplementary figure #4](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s4.png).

## Section 3: Evolutionary analysis

We applied a robust methodology to quantify the dN/dS ratios (non-synonymous to synonymous substitutions) for the annotated genes in the 44 mining MAGs, using the following approach summarised below:

1. Build the reference index with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2):
```
#we index only the binned contigs derived from the S15 metagenome:
bwa-mem2 index S15_bins.fasta
```
Next, we used [alnsl](https://github.com/digenoma-lab/alnsl), an in-house [Nextflow](https://www.nextflow.io/) pipeline which runs ```bwa-mem2``` in order to align the reference genome to the raw reads and [elPrep](https://github.com/ExaScience/elprep) which refines the alignment (e.g., mark duplicated reads), and sorting it to a ```bam``` file (via [samtools](https://github.com/samtools/samtools)). Example execution:
```
#prepare reference file for elPrep:
elprep fasta-to-elfasta S15_bins.fasta S15_bins.fasta.elfasta

#run the pipeline (check both config and environment settings):
nextflow -bg run alnsl/main.nf --csv reads.csv -c alnsl/nextflow.config -profile uoh -params-file alnsl/aln-params.yml
```
The above command gives as the output, among others, a file named ```S15_bins.bam```.

2. Call variants using [inStrain](https://github.com/MrOlm/instrain):
```
#we use the profile module with custom parameters:
inStrain profile --min_mapq 20 --min_read_ani 0.95 S15_bins.bam S15_bins.fasta -o is_S15_profile
```
The variants are contained in the file ```SNVs.tsv```, which is needed to be post-processed to extract and filter relevant information.

3. Install and run [dNdScv](https://github.com/im3sanger/dndscv) to calculate the dN/dS ratios:
```
#first, we parse the GFF file in a tab-delimited CDS table with gene coordinates:
python3 parse_cds.py genes.gff > cds_table.txt

#we now call buildref function (very slow: ~2 hrs in building object):
#here we need the previous .fasta file
buildref(cdsfile = "dndscv/buildref/cds_table.txt",
         genomefile = "dndscv/buildref/S15_bins.fa",
         outfile = "dndscv/buildref/S15_refcds.rda")

#then, we use the SNVs.tsv file to extract the list of mutations (S15_mutations.txt):

#run dNdScv (runtime: ~8 min, high memory usage) and save object:
#this object is later used to calculate dN/dS ratios for each annotated gene or CDS (post-analysis)
snp_dndsout <- dndscv(snp_mutations, refdb = "dndscv/buildref/S15_refcds.rda",
                      max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,
                      cv = NULL, outmats = TRUE)
save(snp_dndsout, file = "dndscv/run/S15_dndsout.rda")

#output objects:
buildref: S15_refcds.rda (118.8 Mb)
dndscv: S15_dndsout.rda (566.9 Mb)
```

See [main figure #3](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_3.png),
[main figure #4](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_4.png).

## Supplementary material

See [supplementary tables](https://docs.google.com/spreadsheets/d/1cci69qkc_zJ21pChGYJtFjSIESJAO0ri/edit?gid=233488985#gid=233488985).

## Directories structure

```
|- data/
| |- dndscv/ (input data for dNdScv package)
| | |- cds_table.txt
| | |- parse_cds.py
| | |- parse_snpeff.py
| | |- S15_mutations.txt
| |- MAG_SMAG/ (MAG tree, taxonomy and coverage data)
| | |- mags_rel_abundances.txt
| | |- mags_S15_r214.txt
| | |- smag_data.tsv
| | |- smag_filtered_bins.tsv
| | |- gtdbtk.DASTool-S15.tree
| |- snps/ (variant calling results and SNPs analysis)
|- figures/
| |- fig_1.png
| |- fig_2.png
| |- fig_3.png
| |- fig_4.png
| |- fig_s1.png
| |- fig_s2.png
| |- fig_s3.png
| |- fig_s4.png
|- R/
| |- paper_filter_metadata.Rmd
| |- paper.Rmd,RData
| |- SNPs.Rmd,html
| |- SNPsv2.Rmd,html
```

## Citation

Moises A. Rojas, Gladis Serrano, Jorge Torres, Jaime Ortega, Gabriel Gálvez, Emilio Vilches, Valentina Parra, Angélica Reyes-Jara, Vinicius Maracaja-Coutinho, Lorena Pizarro, Mauricio Latorre, Alex Di Genova; *"Profiling extremophile bacterial communities recovered from a mining tailing against soil ecosystems through comparative genome-resolved metagenomics and evolutionary analysis"*, bioRxiv, Aug 2024, doi: [10.1101/2024.08.28.610100](https://www.biorxiv.org/content/10.1101/2024.08.28.610100v1).

* Dec 2024: Environmental Microbiome (in revision).

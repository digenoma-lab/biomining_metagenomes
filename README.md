# biomining_metagenomes

Repository that holds data, scripts and figures regarding to Biomining MAGs-SMAG article: **"Profiling extremophile bacterial communities recovered from a mining tailing against soil ecosystems through comparative metagenomics"**.

## Section 1: Taxonomy, QC and abundances

We classified a set of 44 bacterial MAGs belonging to the ```S15``` sample.

See [main figure #1](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_1.png).

See [main figure #2](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_2.png).

See also [supplementary figure #1](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s1.png).

## Section 2: Functional annotation

We analyze our biomining-model ecosystem across conventional ecosystems to compare their capacities over selected pathways (copper, iron, and sulfur). Conventional ecosystems are derived from the [SMAG catalog](https://microbma.github.io/project/SMAG.html) (bacterial subset).

See [main figure #3](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_3.png).

See also [supplementary figure #2](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s2.png).

## Section 3: Evolutionary analysis

1. Build the reference index with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2):
```
#we index only the binned contigs for our S15 sample
bwa-mem2 index S15_bins.fasta
```
Next, we used [alnsl](https://github.com/digenoma-lab/alnsl), a Nextflow pipeline which runs ```bwa-mem2``` to align the reference genome to reads and [elPrep](https://github.com/ExaScience/elprep) to analyze the alignment, refine and sort the resulting ```S15_bins.bam``` file (via [samtools](https://github.com/samtools/samtools))
```
#prepare reference file for elPrep
elprep fasta-to-elfasta ref.fa ref.fa.elfasta
#then, run the pipeline (check both config and environment)
nextflow -bg run alnsl/main.nf --csv reads.csv -c alnsl/nextflow.config -profile uoh -params-file alnsl/aln-params.yml
```
Note: working directory in the Kütral cluster: ```/mnt/beegfs/home/mrojas/DiGenomaLab/Systemix/analysis/snps```
2. Call variants in uncompressed ```VCF``` format with [bcftools](https://samtools.github.io/bcftools/):
```
#mpileup algorithm
bcftools mpileup -Ov -f S15_bins.fasta S15_bins.bam | bcftools call --ploidy 1 -mv -Ov -o S15.vcf
```
Now, we have our ```S15.vcf``` file.
3. Install, configure and run [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/) to annotate variants:
```
#go to dir for preparing files
cd /mnt/beegfs/home/mrojas/DiGenomaLab/Systemix/analysis/snps/snpeff
#copy gff annotations from Prokka
find ../../mags-all-adg/mag-all-samples/Prokka/SPAdes/*Refined-S15*/ -name '*.gff' | xargs cp -t gff/
#link the annotations in the snpEff dir
cd
cd snpEff/data/S15_bins
ln -s ../../../DiGenomaLab/Systemix/analysis/snps/snpeff/gff/genes.gff genes.gff
# Link the fasta reference genome (S15)
cd ../genomes
ln -s ../../../DiGenomaLab/Systemix/analysis/snps/S15_bins.fasta S15_bins.fa

#run SnpEff step #1: build database (takes around 13 min)
#it generates two binary files for our sample: sequence.bin and snpEffectPredictor.bin inside data/S15_bins subdir
java -Xmx16g -jar snpEff.jar build -config snpEff.config -dataDir data/ -noCheckCds -noCheckProtein -gff3 -v S15_bins &> data/S15_bins.build.log

#run SnpEff step #2: annotate (~3 min)
java -Xmx16g -jar snpEff.jar -stats S15.ann.html S15_bins S15.vcf > S15.ann.vcf
```
4. Install and use the R package [dndscv](https://github.com/im3sanger/dndscv) to analyze the annotations:
```
#first, we parse the predicted .gff file in a tab-delimited CDS table with gene coordinates
python3 parse_cds.py genes.gff > cds_table.txt

#call buildref (very slow: ~2 hrs in building object)
#here we need the previously merged/created .fasta file
buildref(cdsfile = "dndscv/buildref/cds_table.txt",
         genomefile = "dndscv/buildref/S15_bins.fa",
         outfile = "dndscv/buildref/S15_refcds.rda")

#then, we use the annotated S15.ann.vcf file to extract mutations
python3 parse_snpeff.py S15.ann.vcf > S15_mutations.txt

#run dndscv (runtime ~8 min, high memory usage) and save object
#this object is later used to calculate dN/dS ratios for each annotated gene or CDS (post-analysis)
snp_dndsout <- dndscv(snp_mutations, refdb = "dndscv/buildref/S15_refcds.rda",
                      max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf,
                      cv = NULL, outmats = TRUE)
save(snp_dndsout, file = "dndscv/run/S15_dndsout.rda")

#resulting objects (stored in the path: /mnt/beegfs/home/mrojas/snpEff/dndscv):
buildref: S15_refcds.rda (118.8 Mb)
output: S15_dndsout.rda (545.9 Mb)
```

See [main figure #4](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_4.png).

See [main figure #5](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_5.png).

See also [supplementary figure #3](https://github.com/digenoma-lab/biomining_metagenomes/blob/main/figures/fig_s3.png).

## Supplementary material

See [supplementary tables](https://docs.google.com/spreadsheets/d/1cci69qkc_zJ21pChGYJtFjSIESJAO0ri/edit?gid=233488985#gid=233488985).

## Directories structure

```
|- data/
| |- dndscv/
| | |- cds_table.txt
| | |- parse_cds.py
| | |- parse_snpeff.py
| | |- S15_DASTool_contig2bin.tsv
| | |- S15_mutations.txt
| |- MAG_SMAG/
| | |- mags_metadata.tsv
| | |- mags_mining_S15.tsv
| | |- smag_filtered_bins.tsv
| |- taxonomy/
| | |- gtdbtk.DASTool-S15.tree
|- figures/
| |- fig_1.png
| |- fig_2.png
| |- fig_3.png
| |- fig_4.png
| |- fig_5.png
| |- fig_s1.png
| |- fig_s2.png
| |- fig_s3.png
|- R/
| |- paper_filter_metadata.Rmd
| |- paper.Rmd
| |- paper.RData
```

## Citation

Moises A. Rojas, Gladis Serrano, Jorge Torres, Jaime Ortega, Gabriel Gálvez, Emilio Vilches, Valentina Parra, Angélica Reyes-Jara, Vinicius Maracaja-Coutinho, Lorena Pizarro, Mauricio Latorre, Alex Di Genova; "Profiling extremophile bacterial communities recovered from a mining tailing against soil ecosystems through comparative genome-resolved metagenomics and evolutionary analysis", bioRxiv, 2024, doi: [10.1101/2024.08.28.610100](https://www.biorxiv.org/content/10.1101/2024.08.28.610100v1).

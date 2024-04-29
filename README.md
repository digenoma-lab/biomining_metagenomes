# biomining_metagenomes

Repository that holds data, scripts and figures regarding to Biomining MAGs-SMAG article: **"Profiling extremophile bacterial communities recovered from a mining tailing against soil ecosystems through comparative metagenomics"**.

## Section 1: Taxonomy

Done.

## Section 2: Metabolism

Done.

## Section 3: Evolution

Working directory in the Kütral cluster: ```/mnt/beegfs/home/mrojas/DiGenomaLab/Systemix/analysis/snps```

**(Current) SNP calling workflow** (```bowtie2/``` subdir)

1. Build the reference index with [bowtie2](https://github.com/BenLangmead/bowtie2):
```
#we use only the contigs that goes to MAG bins for our S15 sample
bowtie2-build S15_bins.fasta S15_build
#then, align genome to reads
bowtie2 -x S15_build -1 S15.R1.fq.gz -2 S15.R2.fq.gz --no-unal -S S15_align.sam
```
2. Convert SAM file to BAM and sort it using [samtools](https://github.com/samtools/samtools):
```
samtools view --threads 8 -bS S15_align.sam -o S15_align.bam
samtools sort --threads 8 S15_align.bam -o S15_sorted.bam
```
3. Generate variant calls in uncompressed ```VCF``` format with [bcftools](https://samtools.github.io/bcftools/):
```
#mpileup algorithm
bcftools mpileup -Ov -f S15_bins.fasta S15_sorted.bam | bcftools call -mv -Ov -o S15.vcf
```
Now, we have our VCF file (```S15.vcf```).
4. Install, configure and run [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/) to annotate variants
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
#run SnpEff step #1: build database (~13 min)
#it generates two binary files for our sample: sequence.bin and snpEffectPredictor.bin inside data/S15_bins subdir
java -Xmx16g -jar snpEff.jar build -config snpEff.config -dataDir data/ -noCheckCds -noCheckProtein -gff3 -v S15_bins &> data/S15_bins.build.log
#run SnpEff step #2: annotate (~3 min)
java -Xmx16g -jar snpEff.jar -stats S15.ann.html S15_bins S15.vcf > S15.ann.vcf
```
5. Install and use the R package [dndscv](https://github.com/im3sanger/dndscv) to analyze the annotations:
> Resulting objects are stored in the path: ```/mnt/beegfs/home/mrojas/snpEff/dndscv```:
```
buildref: S15_refcds.rda
output: S15_dndsout.rda
```


* **To-do (29/04):** Change caller tool to ```bwa-mem2```, which includes quality of SNPs (base quality score recalibration, BQSR).


## Main directories

```
R/
parsers/
checkm_bins.R
get_gtdb_taxa.R
paper_filter_metadata.Rmd
paper.Rmd
paper.RData
```

```
data/
MAG_vs_SMAG/
mags_biomining.tsv
mags_metadata_bac.tsv
mags_metadata.tsv
smag_data.tsv
smag_filtered_bins.tsv
SMAG_mag20177_refined.tree
QC/
checkm_parsed_bins.tsv
checkm_qa.tsv
taxonomy/
bins_refined_list.tsv
gtdbtk_classification.tsv
gtdbtk_summary.tsv
gtdbtk.DASTool-S15.tree
```

```
figures/
final figures in pdf format
```

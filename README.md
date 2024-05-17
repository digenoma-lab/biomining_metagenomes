# biomining_metagenomes

Repository that holds data, scripts and figures regarding to Biomining MAGs-SMAG article: **"Profiling extremophile bacterial communities recovered from a mining tailing against soil ecosystems through comparative metagenomics"**.

## Section 1: Taxonomy

## Section 2: Metabolism

## Section 3: Evolution

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
Now, we have our ```VCF``` file (```S15.vcf```).
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
| | |- mags_biomining.tsv
| | |- mags_metadata.tsv
| | |- mags_metadata_bac.tsv
| | |- smag_data.tsv
| | |- smag_filtered_bins.tsv
| | |- SMAG_mag20177_refined.tree
| |- QC/
| | |- checkm_parsed_bins.tsv
| | |- checkm_qa.tsv
| |- taxonomy/
| | |- bins_refined_list.tsv
| | |- gtdbtk_classification.tsv
| | |- gtdbtk_summary.tsv
| | |- gtdbtk.DASTool-S15.tree
|- figures/
| |- fig_1.pdf
| |- fig_2.pdf
| |- fig_3.png
| |- fig_s1.png
| |- fig_s2.png
| |- fig_s3.png
| |- fig_s3c.png
|- R/
| |- paper_filter_metadata.Rmd
| |- paper.Rmd
| |- paper.RData
```

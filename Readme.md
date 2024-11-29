# Snakemake pipeline for variant calling in bacterial wgs sequences

## Overview
This repository contains a comprehensive bioinformatics pipeline for variant calling in bacterial wgs sequences. The pipeline includes quality control, alignment, variant calling, and annotation steps, and is designed to work with data from Illumina sequencing platforms.

## Prerequisites
Ensure the following tools are installed:
- **Conda** (to manage environments)
- **Snakemake** (for workflow management)
- **Java** (for running tools like snpEff)
- **VCFtools** (for VCF file processing)
- **Other tools** like `bcftools`, `bwa`, `fastqc`, and any additional dependencies defined in the environment YAML files.

## 1. Installation of Snakemake:
   ```bash
   conda install -c bioconda snakemake
   ```

## 2. Sample Handling
**Download the samples**
To download the samples, use the provided getdata.sh script, which automates the download process. This script will pull the necessary FASTQ files based on the sample identifiers in samples.tsv.
```
bash getdata.sh
```
**Sample file for Snakemake**
The pipeline reads sample information from a samples.tsv file. This file should include a list of all sample names and their corresponding file paths. The format should be as follows:
```
sample_name    file_1_path    file_2_path
sample1        data/sample1_1.fastq.gz     data/sample1_2.fastq.gz
sample2        data/sample2_1.fastq.gz     data/sample2_2.fastq.gz
sample3        data/sample3_1.fastq.gz     data/sample3_2.fastq.gz
```

## 3. Preparing the custom snpEff database
**Download the gff file:**
```
curl -L -o data/ref_genome/genes.gff ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.gff.gz
```

**Configure snpEff:**
**Modify the snpEff configuration file (snpEff.config) to include your custom database.**
```
#-------------
# Custom Databases
ecoli_custom.genome=Ecolichromosome
#------------
```
**Copy files into snpEff directory:**
Copy reference fasta and gff into snpEff's path and change file names to sequence.fa and genes.gff. 

**Build snpEff custom database:**
```
snpEff build -gff3 ecoli_custom

```

## 4. Running the Pipeline
**Standard run:**
   ```bash
   snakemake -p --use-conda
   ```

**Optional: Run specific targets** (e.g., alignments, variant calling):
   ```bash
   snakemake -n  # Dry run to see what will be executed
   snakemake <target>  # Execute specific target
   ```

## 5. Outputs
- **Aligned BAM files**: `results/sample1_aligned_sorted.bam`
- **Variant VCF files**: `results/sample1_variants.vcf`
- **Merged VCF file**: `results/cohort.vcf.gz`
- **Annotated VCF**: `results/cohort_annotated.vcf`
- **QC report**: `results/multiqc_report.html`

## 6. Troubleshooting
- Ensure all required software is installed and added to your PATH.
- If you encounter errors related to missing dependencies, activate the relevant conda environment and run the command outside Snakemake.

---

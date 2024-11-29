import pandas as pd

# Specify the configfile
configfile: "config.yaml"

# Load the samples from the TSV file specified in the config.yaml
samples_file = config["samples_file"]
samples= pd.read_csv(samples_file, sep="\t").set_index("sample_id")

#print(samples) # check if samples are correct

# Rule all
rule all:
    input:
        "results/cohort.vcf.gz",
        "results/cohort_annotated.vcf"

# Rule for FastQC
rule fastqc:
    input:
        forward=lambda wildcards: samples.loc[wildcards.sample, "fastq_forward"],
        rev=lambda wildcards: samples.loc[wildcards.sample, "fastq_reverse"],
    output:
        forward_qc="results/{sample}_1_fastqc.html",  # FastQC report for forward read
        reverse_qc="results/{sample}_2_fastqc.html"   # FastQC report for reverse read
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        fastqc {input.forward} --outdir=results
        fastqc {input.rev} --outdir=results
        """

# Rule for MultiQC
rule multiqc:
    input:
        "results/fastqc.html"
    output:
        "results/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        "multiqc -o results {input}"

# Rule for bwa
rule bwa:
    input:
        forward=lambda wildcards: samples.loc[wildcards.sample, "fastq_forward"],
        rev=lambda wildcards: samples.loc[wildcards.sample, "fastq_reverse"]
    output:
        "results/{sample}_aligned.bam"
    conda:
        "bwa_env"
    params:
        ref=config["ref_file"]
    shell:
        """
        bwa mem {params.ref} {input.forward} {input.rev} | samtools view -bS - > {output}
        """


# Rule for Sorting BAM
rule sort_bam:
    input:
        "results/{sample}_aligned.bam"
    output:
        "results/{sample}_aligned_sorted.bam"
    conda:
        "bwa_env"
    shell:
        "samtools sort -o {output} {input}"

# # Rule for Indexing BAM
rule index_bam:
    input:
        "results/{sample}_aligned_sorted.bam"
    output:
        "results/{sample}_aligned_sorted.bam.bai"
    conda:
        "bwa_env"
    shell:
        "samtools index {input}"

# # Rule for Variant Calling
rule variant_calling:
    input:
        bam="results/{sample}_aligned_sorted.bam",
        bai="results/{sample}_aligned_sorted.bam.bai"
    output:
        "results/{sample}_variants.vcf"
    conda:
        "envs/bcftools.yaml"
    params:
        ref=config["ref_file"]
    shell:
        "bcftools mpileup -Ou -f {params.ref} {input.bam} | "
        "bcftools call -mv -Oz -o {output}"


# Compress vcf
rule compress_vcf:
    input:
        "results/{sample}_variants.vcf"
    output:
        "results/{sample}_variants.vcf.gz"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"

# Index vcf
rule index_vcf:
    input:
        "results/{sample}_variants.vcf.gz"
    output:
        "results/{sample}_variants.vcf.gz.csi"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools index {input}"

# Joint genotyping
rule joint_genotyping:
    input:
        vcfs=expand("results/{sample}_variants.vcf.gz", sample=samples.index),
        index=expand("results/{sample}_variants.vcf.gz.csi", sample=samples.index)
    output:
        "results/cohort.vcf.gz"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools merge -o {output} -Oz {input.vcfs}"

# Rule for annotating the merged VCF
rule annotate_merged_vcf:
    input:
        "results/cohort.vcf.gz"
    output:
        "results/cohort_annotated.vcf"
    conda:
        "snpeff_env"
    shell:
        "snpEff ann -v ecoli_custom {input} > {output}"

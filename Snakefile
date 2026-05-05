

import pandas as pd
import os

configfile: "config.yaml"

samples = pd.read_csv(config["samples"]).set_index("sample", drop=False)

def get_srrs_for_gsm(wildcards):
    srrs = samples[samples["gsm"] == wildcards.gsm]["srr"].tolist()
    return expand("data/fastq/{srr}.fastq.gz", srr=srrs)

def get_project_for_gsm(gsm):
    return samples[samples["gsm"] == gsm]["project"].iloc[0]

ALL_GSMS = samples["gsm"].unique().tolist()
ALL_SRRS = samples["srr"].tolist()

rule all:
    input:
        expand("results/fastqc/raw/{srr}_fastqc.html", srr=ALL_SRRS),
        "results/multiqc/raw/multiqc_report.html",
        expand("results/fastqc/trimmed/{gsm}_trimmed_fastqc.html", gsm=ALL_GSMS),
        "results/multiqc/trimmed/multiqc_report.html",
        expand("results/quantifier/{gsm}/miRBase.mrd", gsm=ALL_GSMS),
        "results/counts/merged_counts.csv",
        "results/deseq2/results_PDAC_vs_control.csv",
        "results/deseq2/PCA_plot.pdf",
        "results/deseq2/volcano_plot.pdf"

rule prefetch:
    output:
        sra = "data/sra/{srr}/{srr}.sra"
    params:
        outdir = "data/sra"
    log:
        "logs/prefetch/{srr}.log"
    shell:
        """
        prefetch {wildcards.srr} --output-directory {params.outdir} 2> {log}
        """

rule fasterq_dump:
    input:
        sra = "data/sra/{srr}/{srr}.sra"
    output:
        fastq = temp("data/fastq/{srr}.fastq")
    threads: 4
    log:
        "logs/fasterq_dump/{srr}.log"
    shell:
        """
        fasterq-dump {input.sra} --outdir data/fastq --threads {threads} 2> {log}
        """

rule compress_fastq:
    input:
        fastq = "data/fastq/{srr}.fastq"
    output:
        fastq_gz = "data/fastq/{srr}.fastq.gz"
    shell:
        "gzip -c {input.fastq} > {output.fastq_gz}"

rule fastqc_raw:
    input:
        "data/fastq/{srr}.fastq.gz"
    output:
        html = "results/fastqc/raw/{srr}_fastqc.html",
        zip  = "results/fastqc/raw/{srr}_fastqc.zip"
    threads: 2
    log:
        "logs/fastqc/raw/{srr}.log"
    shell:
        """
        mkdir -p results/fastqc/raw
        fastqc {input} --outdir results/fastqc/raw --threads {threads} 2> {log}
        """

rule multiqc_raw:
    input:
        expand("results/fastqc/raw/{srr}_fastqc.zip", srr=ALL_SRRS)
    output:
        "results/multiqc/raw/multiqc_report.html"
    log:
        "logs/multiqc/raw.log"
    shell:
        """
        multiqc results/fastqc/raw --outdir results/multiqc/raw --force 2> {log}
        """

rule merge_technical_replicates:
    input:
        get_srrs_for_gsm
    output:
        merged = "data/merged/{gsm}.fastq.gz"
    log:
        "logs/merge/{gsm}.log"
    shell:
        """
        cat {input} > {output.merged} 2> {log}
        """


rule trim_galore:
    input:
        "data/merged/{gsm}.fastq.gz"
    output:
        trimmed = "data/trimmed/{gsm}_trimmed.fq.gz"
    params:
        outdir  = "data/trimmed",
        length  = config.get("min_length", 17),
        quality = config.get("min_quality", 20)
    threads: 4
    log:
        "logs/trim_galore/{gsm}.log"
    shell:
        """
        mkdir -p {params.outdir}
        trim_galore \
            --illumina \
            --length {params.length} \
            --quality {params.quality} \
            --cores {threads} \
            --gzip \
            --output_dir {params.outdir} \
            {input} \
            2> {log}
        """

rule fastqc_trimmed:
    input:
        "data/trimmed/{gsm}_trimmed.fq.gz"
    output:
        html = "results/fastqc/trimmed/{gsm}_trimmed_fastqc.html",
        zip  = "results/fastqc/trimmed/{gsm}_trimmed_fastqc.zip"
    threads: 2
    log:
        "logs/fastqc/trimmed/{gsm}.log"
    shell:
        """
        mkdir -p results/fastqc/trimmed
        fastqc {input} --outdir results/fastqc/trimmed --threads {threads} 2> {log}
        """

rule multiqc_trimmed:
    input:
        expand("results/fastqc/trimmed/{gsm}_trimmed_fastqc.zip", gsm=ALL_GSMS)
    output:
        "results/multiqc/trimmed/multiqc_report.html"
    log:
        "logs/multiqc/trimmed.log"
    shell:
        """
        multiqc results/fastqc/trimmed --outdir results/multiqc/trimmed --force 2> {log}
        """

rule collapse_reads:
    input:
        "data/trimmed/{gsm}_trimmed.fq.gz"
    output:
        collapsed = "data/collapsed/{gsm}_collapsed.fa"
    params:
        min_length = config.get("min_length", 17)
    log:
        "logs/collapse/{gsm}.log"
    shell:
        """
        mkdir -p data/collapsed
        zcat {input} > data/trimmed/{wildcards.gsm}_tmp.fq
        /home/user/miniconda3/envs/mirdeep2_env/bin/mapper.pl \
            data/trimmed/{wildcards.gsm}_tmp.fq \
            -e -h -i -j \
            -l {params.min_length} \
            -m \
            -s {output.collapsed} \
            2> {log}
        rm -f data/trimmed/{wildcards.gsm}_tmp.fq
        """

rule quantifier:
    input:
        collapsed    = "data/collapsed/{gsm}_collapsed.fa",
        mature       = config["mirbase_mature_hsa"],
        hairpin      = config["mirbase_hairpin_hsa"]
    output:
        mrd = "results/quantifier/{gsm}/miRBase.mrd"
    params:
        outdir = "results/quantifier/{gsm}"
    log:
        "logs/quantifier/{gsm}.log"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        /home/user/miniconda3/envs/mirdeep2_env/bin/quantifier.pl \
            -p ../../../{input.hairpin} \
            -m ../../../{input.mature} \
            -r ../../../{input.collapsed} \
            -t hsa \
            -y now \
            -d \
            2> ../../../{log}
        mv expression_analyses/expression_analyses_now/miRBase.mrd miRBase.mrd
        """

rule merge_counts:
    input:
        expand("results/quantifier/{gsm}/miRBase.mrd", gsm=ALL_GSMS)
    output:
        "results/counts/merged_counts.csv"
    params:
        samples_csv = config["samples"]
    log:
        "logs/merge_counts.log"
    script:
        "scripts/merge_counts.py"

rule deseq2:
    input:
        counts   = "results/counts/merged_counts.csv",
        metadata = config["samples"]
    output:
        results = "results/deseq2/results_PDAC_vs_control.csv",
        pca     = "results/deseq2/PCA_plot.pdf",
        volcano = "results/deseq2/volcano_plot.pdf"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "scripts/deseq2.R"

import json

FILES = json.load(open(config['SAMPLES_JSON']))

import csv
import os


SAMPLES = sorted(FILES.keys())

MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type+ "_" + assay)

CUT_TAG = config["c_t"]
CHIP = config["chip"]

CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
ALL_SAMPLES =  CHIPS + CUT_TAGS

RUNID = config["RUN_ID"]

BAM=expand("{myrun}/dedup/picard/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/dedup/picard/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG= expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
MACS3 = expand("{myrun}/peaks/macs3/{sample}_peaks.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
MACS3_BAMPE = expand("{myrun}/peaks/macs3/BAMPE/{sample}_peaks.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
SIZE=expand("{myrun}/dedup/picard/insert/{sample}_insert.pdf", sample = ALL_SAMPLES,myrun=RUNID)
BED_DEDUP=expand("{myrun}/bed/bedtools/{sample}.bed", sample = ALL_SAMPLES,myrun=RUNID)
PRESEQ=expand("{myrun}/preseq/{sample}_yield.txt", sample = ALL_SAMPLES,myrun=RUNID)

TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(MACS3)
TARGETS.extend(MACS3_BAMPE)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(SIZE)
TARGETS.extend(BED_DEDUP)
TARGETS.extend(PRESEQ)



rule all:
    input: TARGETS

import re

def get_all_R1(wildcards):
    sample = wildcards.sample
    files = FILES[sample.split('_')[0]][sample.split('_')[1]][sample.split('_')[2]]
    r1_files = [f for f in files if re.search(r'_R1_', f)]
    return sorted(r1_files)

def get_all_R2(wildcards):
    sample = wildcards.sample
    files = FILES[sample.split('_')[0]][sample.split('_')[1]][sample.split('_')[2]]
    r2_files = [f for f in files if re.search(r'_R2_', f)]
    return sorted(r2_files)

rule concatenate_lanes:
    input:
        R1 = get_all_R1,
        R2 = get_all_R2
    output:
        R1_merged = temp("{myrun}/merged/{sample}_R1.fastq.gz"),
        R2_merged = temp("{myrun}/merged/{sample}_R2.fastq.gz")
    params:
        dir = "{myrun}/merged/"
    threads: 1
    shell:
        """
        mkdir -p {params.dir}
        
        cat {input.R1} > {output.R1_merged}
        cat {input.R2} > {output.R2_merged}
        """

rule trimming_trimmomatic:
    input:
        R1 = "{myrun}/merged/{sample}_R1.fastq.gz",
        R2 = "{myrun}/merged/{sample}_R2.fastq.gz"
    output:
        Paired1 = temp("{myrun}/trimmed/{sample}_R1_paired.fastq"),
        Paired2 = temp("{myrun}/trimmed/{sample}_R2_paired.fastq")
    params:
        dir = "{myrun}/trimmed/"
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/cutadapt.yaml"
    shell:
        """
        mkdir -p {params.dir}
        cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -o {output.Paired1} -p {output.Paired2} {input.R1} {input.R2}
        """

rule aligning_bowtie2:
    input:
        Paired1= "{myrun}/trimmed/{sample}_R1_paired.fastq",
        Paired2 = "{myrun}/trimmed/{sample}_R2_paired.fastq",
    output:
        sam = temp("{myrun}/mapped/bowtie2/{sample}.sam")
    params:
        index = config.get('index_bt2_hg','') ,
        dir      = "{myrun}/mapped/bowtie2"
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/bowtie2.yml"
    shell:
        """
        mkdir -p {params.dir}
        bowtie2 -x {params.index} -1 {input.Paired1} -2 {input.Paired2} -S {output.sam} -p {threads} --very-sensitive
        """

rule sorted_samtools:
    input:
        sam = "{myrun}/mapped/bowtie2/{sample}.sam"
    output:
        bam = temp("{myrun}/sorted/samtools/{sample}.bam")
    params:
        dir = "{myrun}/sorted/samtools/"
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        samtools sort -o {output.bam} -O bam {input.sam} -@ {threads}
        samtools index {output.bam}
        """

rule dedup_picard:
    input:
        bam =  "{myrun}/sorted/samtools/{sample}.bam"
    output:
        RG = "{myrun}/dedup/picard/{sample}_RG.bam",
        dedup = "{myrun}/dedup/picard/{sample}.bam",
        bai   = "{myrun}/dedup/picard/{sample}.bam.bai",
        metrics = "{myrun}/dedup/picard/{sample}.bam_metrics.txt"
    params:
        dir="{myrun}/dedup/picard/",
        tmp="{myrun}/dedup/picard/tmp"
    resources:
        mem_mb=140000
    threads: config['THREADS'] 
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        picard AddOrReplaceReadGroups I={input.bam} O={output.RG} RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1
        
        picard MarkDuplicates I={output.RG} O={output.dedup} M={output.metrics} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=coordinate REMOVE_DUPLICATES=true TMP_DIR={params.tmp} 

        samtools index {output.dedup} -@ {threads}

        """

rule filter_chr_samtools:
    """
    Filter out Non-canonical chromosomes, Y and X  and mitochondrial
    """
    input:
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
    output:
        filter = "{myrun}/filter/samtools/{sample}.bam",
        bai="{myrun}/filter/samtools/{sample}.bam.bai"
    params:
        dir    = "{myrun}/filter/samtools/"
    resources:
        mem_mb=140000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}

        samtools view  -@ {threads} -b {input.dedup} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > {output.filter} 
        
        samtools index {output.filter} -@ {threads}
        
        """

rule coverage_deeptools:
    input:
        dedup  = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bw="{myrun}/coverage/deeptools/{sample}_RPKM.bw"
    params:
        norm= config.get('norm_method','RPKM'),
        binsize = config.get('binsize','10'),
        smooth = config.get('smooth_length','300')
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/deeptools.yml"
    shell:
        """
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.dedup} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.binsize} --smoothLength {params.smooth} --numberOfProcessors {threads}
        """

rule macs3:
    input:
        treatment = "{myrun}/filter/samtools/{sample}.bam"
    output:
        peaks_narrow = "{myrun}/peaks/macs3/{sample}_peaks.narrowPeak",
        summits = "{myrun}/peaks/macs3/{sample}_summits.bed",
        peaks_xls = "{myrun}/peaks/macs3/{sample}_peaks.xls"
    params:
        outdir = "{myrun}/peaks/macs3/",
        gsize = config.get('genome_size_bp','hs'),
        qvalue = config.get('peaks_qvalue',0.01)
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/macs3.yml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.treatment} --name {wildcards.sample} --outdir {params.outdir} -f BAM --gsize {params.gsize} -q {params.qvalue} --call-summits 2> {params.outdir}/{wildcards.sample}.log
        """

rule MACS3_BAMPE:
    input:
        treatment = "{myrun}/filter/samtools/BAMPE/{sample}.bam"
    output:
        peaks_narrow = "{myrun}/peaks/macs3/BAMPE/{sample}_peaks.narrowPeak",
        summits = "{myrun}/peaks/macs3/BAMPE/{sample}_summits.bed",
        peaks_xls = "{myrun}/peaks/macs3/BAMPE/{sample}_peaks.xls"
    params:
        outdir = "{myrun}/peaks/macs3/BAMPE/",
        gsize = config.get('genome_size_bp','hs'),
        qvalue = config.get('peaks_qvalue',0.01)
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/macs3.yml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.treatment} --name {wildcards.sample} --outdir {params.outdir} -f BAMPE --gsize {params.gsize} -q {params.qvalue} --call-summits --nomodel --keep-dup al 
        """

rule bam_to_bed:
    input: 
        filter = "{myrun}/filter/samtools/{sample}.bam"
    output:
        bed="{myrun}/bed/bedtools/{sample}.bed"
    params:
        dir  = "{myrun}/bed/bedtools/"
    resources:
        mem_mb=64000
    threads: config['THREADS']
    conda:
        "/home/mattia/miniconda3/envs/bedtools.yml"
    shell:
        """

        mkdir -p {params.dir}

        bedtools bamtobed -i {input.filter}  > {output.bed}

        """
rule stat_samtools:
    input:
        dedup  = "{myrun}/filter/samtools/{sample}.bam"
    output:
        flagstat = "{myrun}/dedup/picard/{sample}.flagstat"
    threads: config.get('THREADS',4)
    conda: "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        samtools flagstat {input.dedup} > {output.flagstat}
        """
rule preseq:
    input:
        bam = "{myrun}/filter/samtools/{sample}.bam"  # Fixed input path
    output:
        preseq = "{myrun}/preseq/{sample}_yield.txt"
    params:
        outdir = "{myrun}/preseq/macs3/"  # Should be 0.01 in config
    resources: 
        mem_mb=64000
    threads: config['THREADS']
    log: "{myrun}/peaks/macs3/{sample}.log"
    conda:
        "/home/mattia/miniconda3/envs/preseq.yml"
    shell:
        """
        mkdir -p {params.outdir}

        preseq lc_extrap -B -P -D -o {output.preseq} {input.bam} 
        """

rule insertsize_picard:
    input:
        dedup = "{myrun}/filter/samtools/{sample}.bam"
    output:
        metrics="{myrun}/dedup/picard/insert/{sample}_insert.txt",
        pdf="{myrun}/dedup/picard/insert/{sample}_insert.pdf"
    params:
        dir="{myrun}/dedup/picard/insert",
        tmp="{myrun}/dedup/picard/tmp"
    resources:
        mem_mb=140000
    threads: config['THREADS'] 
    conda:
        "/home/mattia/miniconda3/envs/samtools.yml"
    shell:
        """
        mkdir -p {params.dir}
        
        picard CollectInsertSizeMetrics I={input.dedup} O={output.metrics} H={output.pdf} M=0.5

        """
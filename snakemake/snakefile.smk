configfile: "config_CutTag.yaml"
import json, re

FILES = json.load(open(config['SAMPLES_JSON']))

SAMPLES = sorted(FILES.keys())

MARK_SAMPLES = []
for sample in SAMPLES:
    for sample_type in FILES[sample].keys():
        for assay in FILES[sample][sample_type].keys():
            MARK_SAMPLES.append(sample + "_" + sample_type+ "_" + assay)

CUT_TAG = config.get("c_t", "CutTag")
CHIP = config.get("chip", "ChIP")

CUT_TAGS  = [sample for sample in MARK_SAMPLES if CUT_TAG in sample]
CHIPS = [sample for sample in MARK_SAMPLES if CHIP in sample]
ALL_SAMPLES =  CHIPS + CUT_TAGS

RUNID = config.get('RUN_ID','run')

BAM=expand("{myrun}/dedup/picard/{sample}.bam", sample=ALL_SAMPLES, myrun=RUNID)
ALL_FLAGSTAT = expand("{myrun}/dedup/picard/{sample}.flagstat", sample = ALL_SAMPLES,myrun=RUNID)
ALL_BIGWIG= expand("{myrun}/coverage/deeptools/{sample}_RPKM.bw", sample = ALL_SAMPLES,myrun=RUNID)
MACS3 = expand("{myrun}/peaks/macs3/{sample}_peaks.narrowPeak", sample = ALL_SAMPLES,myrun=RUNID)
SIZE=expand("{myrun}/dedup/picard/insert/{sample}_insert.pdf", sample = ALL_SAMPLES,myrun=RUNID)
BED_DEDUP=expand("{myrun}/bed/bedtools/dedup/{sample}.bed", sample = ALL_SAMPLES,myrun=RUNID)
PRESEQ=expand("{myrun}/preseq/{sample}_yield.txt", sample = ALL_SAMPLES,myrun=RUNID)

TARGETS = []
TARGETS.extend(BAM)
TARGETS.extend(ALL_BIGWIG)
TARGETS.extend(MACS3)
TARGETS.extend(ALL_FLAGSTAT)
TARGETS.extend(SIZE)
TARGETS.extend(BED_DEDUP)
TARGETS.extend(PRESEQ)


rule all:
    input: TARGETS


def _sample_parts(sample):
    # ensure we operate on the base sample name (strip any accidental path)
    base = sample.split('/')[-1]
    parts = base.split('_')
    if len(parts) < 3:
        raise KeyError(f"Unexpected sample name format: '{sample}' (expected 'sample_sampleType_assay')")
    return parts[0], parts[1], parts[2]


def get_R1(wildcards):
    sample = wildcards.sample
    s, t, a = _sample_parts(sample)
    try:
        files = FILES[s][t][a]
    except KeyError:
        raise KeyError(f"Sample '{s}' with type '{t}' and assay '{a}' not found in SAMPLES_JSON")
    for file in files:
        if re.search(r'_R1_|_1.fastq|_1.fq', file):
            return file
    raise ValueError(f"R1 file not found for {sample}")

def get_R2(wildcards):
    sample = wildcards.sample
    s, t, a = _sample_parts(sample)
    try:
        files = FILES[s][t][a]
    except KeyError:
        raise KeyError(f"Sample '{s}' with type '{t}' and assay '{a}' not found in SAMPLES_JSON")
    for file in files:
        if re.search(r'_R2_|_2.fastq|_2.fq', file):
            return file
    return None


rule trimming_trimmomatic:
    input:
        R1 = get_R1,
        R2 = get_R2,
        adapters = config.get('adapters', '')
    output:
        Paired1 = temp("{myrun}/trimmed/{sample}_R1_paired.fastq"),
        Paired2 = temp("{myrun}/trimmed/{sample}_R2_paired.fastq")
    params:
        dir = "{myrun}/trimmed/"
    threads: config.get('THREADS',4)
    conda: "../envs/trimmomatic.yml"
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
    conda: "../envs/bowtie2.yml"
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
    threads: config.get('THREADS',4)
    conda: "../envs/samtools.yml"
    shell:
        """
        mkdir -p {wildcards.myrun if hasattr(wildcards,'myrun') else ''}
        samtools sort -o {output.bam} -O bam {input.sam} -@ {threads}
        samtools index {output.bam}
        """

rule dedup_picard:
    input:
        bam =  "{myrun}/sorted/samtools/{sample}.bam"
    output:
        dedup = "{myrun}/dedup/picard/{sample}.bam",
        bai   = "{myrun}/dedup/picard/{sample}.bam.bai",
        metrics = "{myrun}/dedup/picard/{sample}.bam_metrics.txt"
    threads: config.get('THREADS',4)
    conda: "../envs/samtools.yml"
    shell:
        """
        mkdir -p $(dirname {output.dedup})
        picard MarkDuplicates I={input.bam} O={output.dedup} M={output.metrics} REMOVE_DUPLICATES=true
        samtools index {output.dedup}
        """

rule coverage_deeptools:
    input:
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
    output:
        bw="{myrun}/coverage/deeptools/{sample}_RPKM.bw"
    params:
        norm= config.get('norm_method','RPKM'),
        binsize = config.get('binsize','10'),
        smooth = config.get('smooth_length','300')
    threads: config.get('THREADS',4)
    conda: "../envs/deeptools.yml"
    shell:
        """
        mkdir -p $(dirname {output.bw})
        bamCoverage -b {input.dedup} --outFileName {output.bw} --normalizeUsing {params.norm} --binSize {params.binsize} --smoothLength {params.smooth} --numberOfProcessors {threads}
        """

rule macs3:
    input:
        treatment = "{myrun}/dedup/picard/{sample}.bam"
    output:
        peaks_narrow = "{myrun}/peaks/macs3/{sample}_peaks.narrowPeak",
        summits = "{myrun}/peaks/macs3/{sample}_summits.bed",
        peaks_xls = "{myrun}/peaks/macs3/{sample}_peaks.xls"
    params:
        outdir = "{myrun}/peaks/macs3/",
        gsize = config.get('genome_size_bp','hs'),
        qvalue = config.get('peaks_qvalue',0.01)
    threads: config.get('THREADS',4)
    conda: "../envs/macs3.yml"
    shell:
        """
        mkdir -p {params.outdir}
        macs3 callpeak -t {input.treatment} --name {wildcards.sample} --outdir {params.outdir} -f BAM --gsize {params.gsize} -q {params.qvalue} --call-summits 2> {params.outdir}/{wildcards.sample}.log
        """

rule stat_samtools:
    input:
        dedup  = "{myrun}/dedup/picard/{sample}.bam"
    output:
        flagstat = "{myrun}/dedup/picard/{sample}.flagstat"
    threads: config.get('THREADS',4)
    conda: "../envs/samtools.yml"
    shell:
        """
        samtools flagstat {input.dedup} > {output.flagstat}
        """

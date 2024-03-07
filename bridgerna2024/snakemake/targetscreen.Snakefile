from glob import glob
from os.path import basename, join
from bridgerna2024 import *
from bridgerna2024 import targetscreen

WD = config['workdir']
AMPLICON_DIR = join(WD, config['amplicon_dir'])
FASTQ_DIR = join(WD, config['fastq_dir'])
CLEAN_DIR = join(WD, config['clean_dir'])
AMPLICON_ALIGN_DIR = join(WD, config['amplicon_align_dir'])
BARCODE_DIR = join(WD, config['barcode_dir'])
OLIGO_COVERAGE_DIR = join(WD, config['oligo_coverage_dir'])
VISUALIZE_DIR = join(WD, config['visualize_dir'])

OLIGO_INFO = join(WD, config['oligo_info'])

SAMPLES = list(set([basename(f).split('.')[0] for f in glob(join(FASTQ_DIR, '*fq.gz'))]))

rule all:
    input:
        join(VISUALIZE_DIR, 'cpm_biorep_comparison.pdf'),
        expand(join(CLEAN_DIR, '{sample}.R1_fastqc.html'), sample=SAMPLES)

rule bbduk:
    input:
        fq1=join(FASTQ_DIR, '{sample}.R1.fq.gz'),
        fq2=join(FASTQ_DIR, '{sample}.R2.fq.gz')
    output:
        fq1=join(CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(CLEAN_DIR, '{sample}.R2.fq.gz')
    shell:
        """
        bbduk.sh in1={input.fq1} in2={input.fq2} out1={output.fq1} out2={output.fq2} ref={BBDUK_ADAPTERS_PATH} ktrim=r mink=11 hdist=1 k=23 tpe tbo
        """

rule fastqc:
    input:
        fq1=join(CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(CLEAN_DIR, '{sample}.R2.fq.gz')
    output:
        html1=join(CLEAN_DIR, '{sample}.R1_fastqc.html'),
        html2=join(CLEAN_DIR, '{sample}.R2_fastqc.html')
    threads: 32
    shell:
        """
        fastqc {input.fq1} {input.fq2} --threads {threads}
        """

rule bwa_index:
    input:
        join(AMPLICON_DIR, '{sample}.fasta')
    output:
        join(AMPLICON_DIR, '{sample}.fasta.bwt')
    shell:
        """
        {BWA_PATH} index {input}
        """

rule align_amplicon_bwa:
    input:
        amplicon=join(AMPLICON_DIR, '{sample}.fasta'),
        amplicon_idx=join(AMPLICON_DIR, '{sample}.fasta.bwt'),
        fq1=join(CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(CLEAN_DIR, '{sample}.R2.fq.gz')
    output:
        bam=join(AMPLICON_ALIGN_DIR, '{sample}.bam'),
    threads: 8
    shell:
        """
        {BWA_PATH} mem -t {threads} -T 20 -B 1 -L 500,500 -k 6 {input.amplicon} {input.fq1} {input.fq2} | samtools view -@ {threads} -b > {output.bam};
        """

rule extract_barcodes:
    input:
        bam=join(AMPLICON_ALIGN_DIR, '{sample}.bam'),
        amplicon=join(AMPLICON_DIR, '{sample}.fasta')
    output:
        tsv=join(BARCODE_DIR, '{sample}.barcodes.tsv')
    run:
        targetscreen.extract_barcodes(input.bam, input.amplicon, output.tsv)

rule map_barcodes:
    input:
        tsv=join(BARCODE_DIR, '{sample}.barcodes.tsv'),
    output:
        tsv=join(BARCODE_DIR, '{sample}.barcodes.mapped.tsv')
    run:
        targetscreen.map_barcodes(input.tsv, OLIGO_INFO, output.tsv)

rule oligo_coverage:
    input:
        expand(join(BARCODE_DIR, '{sample}.barcodes.mapped.tsv'), sample=SAMPLES)
    output:
        join(OLIGO_COVERAGE_DIR, 'all_oligo_coverage.tsv')
    params:
        barcode_dir=BARCODE_DIR,
        oligo_info=OLIGO_INFO,
        oligo_coverage_dir=OLIGO_COVERAGE_DIR
    shell:
        """
        Rscript {TARGETSCREEN_OLIGO_COVERAGE_SCRIPT_PATH} {params.barcode_dir} {params.oligo_info} {params.oligo_coverage_dir}
        """

rule visualize:
    input:
        join(OLIGO_COVERAGE_DIR, 'all_oligo_coverage.tsv')
    output:
        join(VISUALIZE_DIR, 'cpm_biorep_comparison.pdf')
    params:
        oligo_info=OLIGO_INFO,
        visualize_dir=VISUALIZE_DIR
    shell:
        """
        Rscript {TARGETSCREEN_VISUALIZE_SCRIPT_PATH} {input} {params.oligo_info} {params.visualize_dir}
        """
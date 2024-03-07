from glob import glob
from os.path import basename, join
from collections import defaultdict
from bridgerna2024 import *
from bridgerna2024 import donorscreen
import sys

WD = config['workdir']

UMI_AMPLICON_DIR = join(WD, config['umi_amplicon_dir'])
UMI_FASTQ_DIR = join(WD, config['umi_fastq_dir'])
UMI_CLEAN_DIR = join(WD, config['umi_clean_dir'])
UMI_ALIGN_DIR = join(WD, config['umi_align_dir'])
UMI_VARIABLE_REGIONS_DIR = join(WD, config['umi_variable_regions_dir'])
UMI_FINAL_MAPPING_DIR = join(WD, config['umi_final_mapping_dir'])

SCREEN_AMPLICON_DIR = join(WD, config['screen_amplicon_dir'])
SCREEN_FASTQ_DIR = join(WD, config['screen_fastq_dir'])
SCREEN_CLEAN_DIR = join(WD, config['screen_clean_dir'])
SCREEN_ALIGN_DIR = join(WD, config['screen_align_dir'])
SCREEN_UMI_DIR = join(WD, config['screen_umi_dir'])
SCREEN_OLIGO_COUNTS_DIR = join(WD, config['screen_oligo_counts_dir'])
VISUALIZE_DIR = join(WD, config['visualize_dir'])

UMI_SAMPLES = list(set([basename(f).split('.')[0] for f in glob(join(UMI_FASTQ_DIR, '*fq.gz'))]))
SCREEN_SAMPLES = list(set([basename(f).split('.')[0] for f in glob(join(SCREEN_FASTQ_DIR, '*fq.gz'))]))

rule all:
    input:
        join(VISUALIZE_DIR, 'cpm_biorep_comparison.pdf'),
        expand(join(UMI_CLEAN_DIR, '{sample}.R1_fastqc.html'), sample=UMI_SAMPLES),
        expand(join(SCREEN_CLEAN_DIR, '{sample}.R1_fastqc.html'), sample=SCREEN_SAMPLES)

rule umi_bbduk_fastq:
    input:
        fq1=join(UMI_FASTQ_DIR, '{sample}.R1.fq.gz'),
        fq2=join(UMI_FASTQ_DIR, '{sample}.R2.fq.gz')
    output:
        fq1=join(UMI_CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(UMI_CLEAN_DIR, '{sample}.R2.fq.gz')
    shell:
        """
        bbduk.sh in1={input.fq1} in2={input.fq2} out1={output.fq1} out2={output.fq2} ref={BBDUK_ADAPTERS_PATH} ktrim=r mink=11 hdist=1 k=23 tpe tbo
        """

rule umi_fastqc_fastq:
    input:
        fq1=join(UMI_CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(UMI_CLEAN_DIR, '{sample}.R2.fq.gz')
    output:
        html1=join(UMI_CLEAN_DIR, '{sample}.R1_fastqc.html'),
        html2=join(UMI_CLEAN_DIR, '{sample}.R2_fastqc.html')
    threads: 32
    shell:
        """
        fastqc {input.fq1} {input.fq2} --threads {threads}
        """

rule umi_amplicon_bwa_index:
    input:
        join(UMI_AMPLICON_DIR, '{sample}.fasta')
    output:
        join(UMI_AMPLICON_DIR, '{sample}.fasta.bwt')
    shell:
        """
        {BWA_PATH} index {input}
        """

rule umi_align_bwa_mem:
    input:
        fq1=join(UMI_CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(UMI_CLEAN_DIR, '{sample}.R2.fq.gz'),
        umi_amp_fasta=join(UMI_AMPLICON_DIR, '{sample}.fasta'),
        umi_amp_idx=join(UMI_AMPLICON_DIR, '{sample}.fasta.bwt')
    output:
        bam=join(UMI_ALIGN_DIR, '{sample}.bam')
    threads: 32
    shell:
        """
        {BWA_PATH} mem -t {threads} -T 20 -B 1 -L 500,500 -k 6 {input.umi_amp_fasta} {input.fq1} {input.fq2} | samtools view -@ {threads} -b > {output.bam}
        """

rule umi_extract_variable_regions:
    input:
        bam=join(UMI_ALIGN_DIR, '{sample}.bam'),
        umi_amp_fasta = join(UMI_AMPLICON_DIR, '{sample}.fasta'),
        umi_amp_bed = join(UMI_AMPLICON_DIR, '{sample}.bed')
    output:
        umis=join(UMI_VARIABLE_REGIONS_DIR, '{sample}.variable_regions.tsv')
    params:
        outdir = UMI_VARIABLE_REGIONS_DIR
    run:
        donorscreen.umi_align_extract_variable_regions(input.bam, input.umi_amp_fasta, input.umi_amp_bed, params.outdir)

rule umi_map_oligos:
    input:
        umis=join(UMI_VARIABLE_REGIONS_DIR, '{sample}.variable_regions.tsv')
    output:
        umis=join(UMI_VARIABLE_REGIONS_DIR, '{sample}.variable_regions.mapped_oligos.tsv')
    params:
        oligos = join(WD, config['oligo_info'])
    run:
        donorscreen.map_oligos(input.umis, params.oligos, output.umis)

rule umi_final_mapping:
    input:
        umis=[join(UMI_VARIABLE_REGIONS_DIR, '{sample}.variable_regions.mapped_oligos.tsv'.format(sample=samp)) for samp in UMI_SAMPLES],
        oligos=join(WD, config['oligo_info'])
    output:
        final=join(UMI_FINAL_MAPPING_DIR, 'mapped_oligos.final.tsv')
    run:
        donorscreen.final_umi_mapping(input.umis, input.oligos, output.final)

rule screen_bbduk_fastq:
    input:
        fq1=join(SCREEN_FASTQ_DIR, '{sample}.R1.fq.gz'),
        fq2=join(SCREEN_FASTQ_DIR, '{sample}.R2.fq.gz')
    output:
        fq1=join(SCREEN_CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(SCREEN_CLEAN_DIR, '{sample}.R2.fq.gz')
    shell:
        """
        bbduk.sh in1={input.fq1} in2={input.fq2} out1={output.fq1} out2={output.fq2} ref={BBDUK_ADAPTERS_PATH} ktrim=r mink=11 hdist=1 k=23 tpe tbo
        """

rule screen_fastqc_fastq:
    input:
        fq1=join(SCREEN_CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(SCREEN_CLEAN_DIR, '{sample}.R2.fq.gz')
    output:
        html1=join(SCREEN_CLEAN_DIR, '{sample}.R1_fastqc.html'),
        html2=join(SCREEN_CLEAN_DIR, '{sample}.R2_fastqc.html')
    threads: 32
    shell:
        """
        fastqc {input.fq1} {input.fq2} --threads {threads}
        """

rule screen_amplicon_bwa_index:
    input:
        join(SCREEN_AMPLICON_DIR, '{sample}.fasta')
    output:
        join(SCREEN_AMPLICON_DIR, '{sample}.fasta.bwt')
    shell:
        """
        {BWA_PATH} index {input}
        """

rule screen_align_amplicon_bwa:
    input:
        amplicon=join(SCREEN_AMPLICON_DIR, '{sample}.fasta'),
        amplicon_idx=join(SCREEN_AMPLICON_DIR, '{sample}.fasta.bwt'),
        fq1=join(SCREEN_CLEAN_DIR, '{sample}.R1.fq.gz'),
        fq2=join(SCREEN_CLEAN_DIR, '{sample}.R2.fq.gz')
    output:
        bam=join(SCREEN_ALIGN_DIR, '{sample}.bam')
    threads: 32
    shell:
        """
        {BWA_PATH} mem -t {threads} -T 20 -B 1 -L 500,500 -k 6 {input.amplicon} {input.fq1} {input.fq2} | samtools view -@ {threads} -b > {output.bam}
        """

rule screen_extract_umis:
    input:
        bam=join(SCREEN_ALIGN_DIR, '{sample}.bam'),
        amplicon=join(SCREEN_AMPLICON_DIR, '{sample}.fasta'),
        amplicon_bed=join(SCREEN_AMPLICON_DIR, '{sample}.bed')
    output:
        tsv=join(SCREEN_UMI_DIR, '{sample}.screen_umis.tsv')
    run:
        donorscreen.screen_extract_umi(input.bam, input.amplicon, input.amplicon_bed, output.tsv)

rule screen_final_oligo_counts:
    input:
        screen_umis=[join(SCREEN_UMI_DIR, '{sample}.screen_umis.tsv'.format(sample=samp)) for samp in SCREEN_SAMPLES],
        umi_map=join(UMI_FINAL_MAPPING_DIR, 'mapped_oligos.final.tsv')
    output:
        counts=join(SCREEN_OLIGO_COUNTS_DIR, 'screen_oligo_counts.tsv')
    run:
        donorscreen.screen_umi_counts(input.screen_umis, input.umi_map, output.counts)

rule visualize:
    input:
        join(SCREEN_OLIGO_COUNTS_DIR, 'screen_oligo_counts.tsv')
    output:
        join(VISUALIZE_DIR, 'cpm_biorep_comparison.pdf')
    params:
        oligo_info=join(WD, config['oligo_info']),
        single_mismatch_sets=join(WD, config['single_mismatch_sets']),
        visualize_dir=VISUALIZE_DIR
    shell:
        """
        Rscript {DONORSCREEN_VISUALIZE_SCRIPT_PATH} {input} {params.oligo_info} {params.single_mismatch_sets} {params.visualize_dir}
        """

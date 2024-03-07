from glob import glob
from os.path import basename, join
from bridgerna2024 import *
from bridgerna2024 import rnaseq

WD = config['workdir']
FASTQ_DIR = join(WD, config['fastq_dir'])
PLASMID_DIR = join(WD, config['plasmid_dir'])
CLEAN_DIR = join(WD, config['clean_dir'])
MERGE_DIR = join(WD, config['merge_dir'])
ALIGN_DIR = join(WD, config['align_dir'])
COVERAGE_DIR = join(WD, config['coverage_dir'])
RESULTS_DIR = join(WD, config['results_dir'])

LIBRARIES = list(set([basename(f).split('.')[0] for f in glob(join(FASTQ_DIR, '*fq.gz'))]))

rule all:
    input:
        join(RESULTS_DIR, 'plasmid_coverage.pdf'),
        expand(join(CLEAN_DIR, '{library}.cleaned.R1_fastqc.html'), library=LIBRARIES)

rule bbduk:
    input:
        fq1=join(FASTQ_DIR, '{library}.R1.fq.gz'),
        fq2=join(FASTQ_DIR, '{library}.R2.fq.gz')
    output:
        fq1=join(CLEAN_DIR, '{library}.cleaned.R1.fq.gz'),
        fq2=join(CLEAN_DIR, '{library}.cleaned.R2.fq.gz')
    shell:
        """
        bbduk.sh in1={input.fq1} in2={input.fq2} out1={output.fq1} out2={output.fq2} ref={BBDUK_ADAPTERS_PATH} ktrim=r mink=11 hdist=1 k=23 tpe tbo
        """

rule bbmerge:
    input:
        fq1=join(CLEAN_DIR, '{library}.cleaned.R1.fq.gz'),
        fq2=join(CLEAN_DIR, '{library}.cleaned.R2.fq.gz')
    output:
        fq1=join(MERGE_DIR, '{library}.unmerged.R1.fq.gz'),
        fq2=join(MERGE_DIR, '{library}.unmerged.R2.fq.gz'),
        fq3=join(MERGE_DIR, '{library}.merged.fq.gz')
    shell:
        """
        bbmerge.sh in1={input.fq1} in2={input.fq2} out={output.fq3} outu1={output.fq1} outu2={output.fq2}
        """

rule bwa_index:
    input:
        fasta=join(PLASMID_DIR, '{library}.fasta')
    output:
        index=join(PLASMID_DIR, '{library}.fasta.bwt')
    shell:
        """
        {BWA_PATH} index {input.fasta}
        """

rule bwa_mem:
    input:
        fq1=join(MERGE_DIR, '{library}.unmerged.R1.fq.gz'),
        fq2=join(MERGE_DIR, '{library}.unmerged.R2.fq.gz'),
        fq3=join(MERGE_DIR, '{library}.merged.fq.gz'),
        fasta=join(PLASMID_DIR, '{library}.fasta'),
        index=join(PLASMID_DIR, '{library}.fasta.bwt')
    output:
        bam1=join(ALIGN_DIR, '{library}.merged.bam'),
        bam2=join(ALIGN_DIR, '{library}.unmerged.bam')
    shell:
        """
        {BWA_PATH} mem -T 20 -t 1 {input.fasta} {input.fq3} | samtools sort > {output.bam1}
        {BWA_PATH} mem -T 20 -t 1 {input.fasta} {input.fq1} {input.fq2} | samtools sort > {output.bam2}
        samtools index {output.bam1}
        samtools index {output.bam2}
        """

rule align_plasmids:
    input:
        all_plasmids=[f for f in glob(join(PLASMID_DIR, '*fasta'))]
    output:
        fasta=join(PLASMID_DIR, 'plasmids_muscle_aligned.fasta')
    shell:
        """
        cat {input.all_plasmids} > tmp.fasta
        muscle -align tmp.fasta -output {output.fasta}; rm tmp.fasta
        """

rule map_alignment_positions:
    input:
        fasta=join(PLASMID_DIR, 'plasmids_muscle_aligned.fasta')
    output:
        tsv=join(PLASMID_DIR, 'plasmids_muscle_aligned.tsv')
    run:
        rnaseq.map_alignment_positions(input.fasta, output.tsv)

rule fastqc:
    input:
        fq1=join(CLEAN_DIR, '{library}.cleaned.R1.fq.gz'),
        fq2=join(CLEAN_DIR, '{library}.cleaned.R2.fq.gz')
    output:
        html1=join(CLEAN_DIR, '{library}.cleaned.R1_fastqc.html'),
        html2=join(CLEAN_DIR, '{library}.cleaned.R2_fastqc.html')
    shell:
        """
        fastqc {input.fq1} {input.fq2} --threads 1
        """

rule calculate_coverage:
    input:
        bam1=join(ALIGN_DIR, '{library}.merged.bam'),
        bam2=join(ALIGN_DIR, '{library}.unmerged.bam'),
        plasmid_fasta=join(PLASMID_DIR, '{library}.fasta'),
        plasmid_bed=join(PLASMID_DIR, '{library}.bed')
    output:
        depth=join(COVERAGE_DIR, '{library}.read_depth.tsv'),
        abundance=join(COVERAGE_DIR, '{library}.read_abundance.tsv')
    run:
        rnaseq.calculate_coverage(
        input.bam1, input.bam2, input.plasmid_fasta, input.plasmid_bed, output.depth, output.abundance
        )

rule visualize:
    input:
        indepth=expand(join(COVERAGE_DIR, '{library}.read_depth.tsv'), library=LIBRARIES),
        intsv=join(PLASMID_DIR, 'plasmids_muscle_aligned.tsv')
    output:
        join(RESULTS_DIR, 'plasmid_coverage.pdf')
    shell:
        """
        Rscript {RNASEQ_VISUALIZE_SCRIPT_PATH} {COVERAGE_DIR} {PLASMID_DIR} {input.intsv} {output}
        """
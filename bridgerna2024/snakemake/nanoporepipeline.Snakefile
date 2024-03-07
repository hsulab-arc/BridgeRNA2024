from glob import glob
from os.path import basename, join
from bridgerna2024 import *
from bridgerna2024 import nanoporepipeline

WD = config['workdir']
FASTQ_DIR = join(WD, config['fastq_dir'])
GENOME_DIR = join(WD, config['genome_dir'])
PLASMID_DIR = join(WD, config['plasmid_dir'])
BRIDGERNA_DIR = join(WD, config['bridgerna_dir'])

CHECK_INDEX_HOPPING_DIR = join(WD, config['check_index_hopping_dir'])
COMBINED_DOWNSAMPLED_FASTQ_DIR = join(WD, config['combined_downsampled_fastq_dir'])
FULL_GENOME_ALIGN_DIR = join(WD, config['full_genome_align_dir'])
FULL_PLASMID_ALIGN_DIR = join(WD, config['full_plasmid_align_dir'])
EXTRACTED_FLANKS_DIR = join(WD, config['extracted_flanks_dir'])
PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR = join(WD, config['plasmid_aligned_extracted_flanks_dir'])
GENOME_ALIGNED_EXTRACTED_FLANKS_DIR = join(WD, config['genome_aligned_extracted_flanks_dir'])
GENOME_ALIGN_NOPLASMID_DIR = join(WD, config['genome_align_noplasmid_dir'])
INSERTION_SITES_DIR = join(WD, config['insertion_sites_dir'])
RAW_STRUCTURAL_VARIANTS_DIR = join(WD, config['raw_structural_variants_dir'])
INSERTION_ANALYSIS_DIR = join(WD, config['insertion_analysis_dir'])
FINAL_STRUCTURAL_VARIANTS_DIR = join(WD, config['final_structural_variants_dir'])

GENOME_FASTA = join(GENOME_DIR, config['genome_fasta'])
SAMPLE_INFO = join(WD, config['sample_info'])

samples_biorep_techrep = list(set(['.'.join(basename(x).split('.')[:-2]) for x in glob(join(FASTQ_DIR, '*.fq.gz'))]))
samples_biorep = list(set(['.'.join(x.split('.')[:-1]) for x in samples_biorep_techrep]))
genomes = set([basename(x).split('.')[0] for x in glob(join(GENOME_DIR, '*.fasta'))])

rule all:
    input:
        expand(join(FASTQ_DIR, '{sample}.read_lengths.tsv'), sample=samples_biorep_techrep),
        join(CHECK_INDEX_HOPPING_DIR, 'combined_bridgerna_guide_counts.tsv'),
        expand(join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.read_lengths.tsv'), sample=samples_biorep),
        join(FINAL_STRUCTURAL_VARIANTS_DIR, 'final_indels.tsv')

rule raw_read_length_distribution:
    input:
        fastq=join(FASTQ_DIR, '{sample}.fq.gz')
    output:
        tsv=join(FASTQ_DIR, '{sample}.read_lengths.tsv')
    run:
        nanoporepipeline.read_length_distribution(input.fastq, output.tsv)

rule bridgerna_alignment:
    input:
        idx=join(BRIDGERNA_DIR, 'bridgerna.mmi'),
        fastq=join(FASTQ_DIR, '{sample}.fq.gz')
    output:
        bam=join(CHECK_INDEX_HOPPING_DIR, '{sample}.bam'),
        bai=join(CHECK_INDEX_HOPPING_DIR, '{sample}.bam.bai')
    threads: 8
    shell:
        """
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.fastq} | samtools view -b -F 4 | samtools sort -@ {threads} > {output.bam};
        samtools index {output.bam}
        """

rule parse_bridgerna_guides:
    input:
        bam=join(CHECK_INDEX_HOPPING_DIR, '{sample}.bam'),
        fasta=join(BRIDGERNA_DIR, 'bridgerna.fasta')
    output:
        guides=join(CHECK_INDEX_HOPPING_DIR, '{sample}.guides.tsv')
    run:
        nanoporepipeline.parse_bridgerna_guides(input.bam, input.fasta, output.guides)

rule count_bridgerna_guides:
    input:
        guides=expand(join(CHECK_INDEX_HOPPING_DIR, '{sample}.guides.tsv'), sample=samples_biorep_techrep)
    output:
        counts=join(CHECK_INDEX_HOPPING_DIR, 'combined_bridgerna_guide_counts.tsv')
    run:
        nanoporepipeline.count_bridgerna_guides(input.guides, output.counts)

rule combine_downsample_fastq:
    output:
        fastq=join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.{biorep}.fq.gz')
    params:
        sample='{sample}',
        biorep='{biorep}',
        sample_info=SAMPLE_INFO
    run:
        infastqs = glob(join(FASTQ_DIR, f'{params.sample}.{params.biorep}.*fq.gz'))
        nanoporepipeline.combine_downsample_fastq(params.sample, infastqs, params.sample_info, output.fastq)

rule combined_downsampled_read_length_distribution:
    input:
        fastq=join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.fq.gz')
    output:
        tsv=join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.read_lengths.tsv')
    run:
        nanoporepipeline.read_length_distribution(input.fastq, output.tsv)

rule minimap2_index_genome:
    input:
        GENOME_FASTA
    output:
        GENOME_FASTA.replace('.fasta', '.mmi')
    shell:
        """
        minimap2 -d {output} {input}
        """

rule minimap2_index_plasmid:
    input:
        join(PLASMID_DIR, '{sample}.{biorep}.fasta')
    output:
        join(PLASMID_DIR, '{sample}.{biorep}.mmi')
    shell:
        """
        minimap2 -d {output} {input}
        """

rule full_genome_alignment:
    input:
        idx=GENOME_FASTA.replace('.fasta', '.mmi'),
        fastq=join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.{biorep}.fq.gz')
    output:
        bam=join(FULL_GENOME_ALIGN_DIR, '{sample}.{biorep}.bam'),
        bai=join(FULL_GENOME_ALIGN_DIR, '{sample}.{biorep}.bam.bai')
    threads: 16
    shell:
        """
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.fastq} | samtools sort -@ {threads} > {output.bam};
        samtools index {output.bam}
        """

rule full_plasmid_alignment:
    input:
        idx=join(PLASMID_DIR, '{sample}.{biorep}.mmi'),
        fastq=join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.{biorep}.fq.gz')
    output:
        bam=join(FULL_PLASMID_ALIGN_DIR, '{sample}.{biorep}.bam'),
        bai=join(FULL_PLASMID_ALIGN_DIR, '{sample}.{biorep}.bam.bai')
    threads: 16
    shell:
        """
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.fastq} | samtools sort -@ {threads} > {output.bam};
        samtools index {output.bam}
        """

rule extract_flanks_from_donor_reads:
    input:
        fastq=join(COMBINED_DOWNSAMPLED_FASTQ_DIR, '{sample}.{biorep}.fq.gz'),
    output:
        left_flank=join(EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.left_flank.fastq'),
        right_flank=join(EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.right_flank.fastq')
    params:
        sample_id='{sample}',
        sample_info=SAMPLE_INFO
    run:
        nanoporepipeline.extract_flanks_from_donor_reads(params.sample_id, input.fastq, params.sample_info, output.left_flank, output.right_flank)

rule align_flanks_to_plasmid:
    input:
        idx=join(PLASMID_DIR, '{sample}.{biorep}.mmi'),
        leftflanks=join(EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.left_flank.fastq'),
        rightflanks=join(EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.right_flank.fastq')
    output:
        leftflanks_bam=join(PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.leftflanks.bam'),
        leftflanks_bai=join(PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.leftflanks.bam.bai'),
        rightflanks_bam=join(PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.rightflanks.bam'),
        rightflanks_bai=join(PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.rightflanks.bam.bai')
    threads: 8
    shell:
        """
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.leftflanks} | samtools sort -@ {threads} > {output.leftflanks_bam};
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.rightflanks} | samtools sort -@ {threads} > {output.rightflanks_bam};
        samtools index {output.leftflanks_bam}
        samtools index {output.rightflanks_bam}
        """

rule align_flanks_to_genome:
    input:
        idx=GENOME_FASTA.replace(".fasta", ".mmi"),
        leftflanks=join(EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.left_flank.fastq'),
        rightflanks=join(EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.right_flank.fastq')
    output:
        leftflanks_bam=join(GENOME_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.leftflanks.bam'),
        leftflanks_bai=join(GENOME_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.leftflanks.bam.bai'),
        rightflanks_bam=join(GENOME_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.rightflanks.bam'),
        rightflanks_bai=join(GENOME_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.rightflanks.bam.bai')
    threads: 8
    shell:
        """
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.leftflanks} | samtools sort -@ {threads} > {output.leftflanks_bam};
        minimap2 --MD -m 16 -ax map-ont -t {threads} {input.idx} {input.rightflanks} | samtools sort -@ {threads} > {output.rightflanks_bam};
        samtools index {output.leftflanks_bam}
        samtools index {output.rightflanks_bam}
        """

rule identify_insertion_sites:
    input:
        genome_leftflanks_bam=join(GENOME_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.leftflanks.bam'),
        genome_rightflanks_bam=join(GENOME_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.rightflanks.bam'),
        plasmid_leftflanks_bam=join(PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.leftflanks.bam'),
        plasmid_rightflanks_bam=join(PLASMID_ALIGNED_EXTRACTED_FLANKS_DIR, '{sample}.{biorep}.rightflanks.bam'),
        genome_fasta=GENOME_FASTA
    output:
        insertion_sites_bam=join(INSERTION_SITES_DIR, '{sample}.{biorep}.insertion_sites.bam'),
        insertion_sites_bai=join(INSERTION_SITES_DIR, '{sample}.{biorep}.insertion_sites.bam.bai'),
        insertion_sites_tsv=join(INSERTION_SITES_DIR, '{sample}.{biorep}.insertion_sites.tsv')
    run:
        nanoporepipeline.identify_insertion_sites(input.genome_leftflanks_bam, input.genome_rightflanks_bam,
        input.plasmid_leftflanks_bam, input.plasmid_rightflanks_bam, input.genome_fasta, output.insertion_sites_bam,
        output.insertion_sites_tsv)

rule merge_insertion_sites:
    input:
        insertion_sites_tsv=join(INSERTION_SITES_DIR, '{sample}.{biorep}.insertion_sites.tsv')
    output:
        merged_insertion_sites_tsv=join(INSERTION_SITES_DIR, '{sample}.{biorep}.insertion_sites.merged.tsv')
    run:
        nanoporepipeline.merge_insertion_sites(input.insertion_sites_tsv, output.merged_insertion_sites_tsv)

rule final_insertion_analysis:
    input:
        merged_insertion_sites_tsv=expand(join(INSERTION_SITES_DIR, '{sample}.insertion_sites.merged.tsv'), sample=samples_biorep)
    params:
        sample_info=SAMPLE_INFO,
        insertion_sites_dir=INSERTION_SITES_DIR,
        outdir=INSERTION_ANALYSIS_DIR
    output:
        join(INSERTION_ANALYSIS_DIR, 'programmed_insertion_sites.tsv'),
        join(INSERTION_ANALYSIS_DIR, 'wt_insertion_sites.tsv')
    shell:
        """
        Rscript {NANOPOREPIPELINE_INSERTION_ANALYSIS_SCRIPT_PATH} {params.insertion_sites_dir} {params.sample_info} {params.outdir}
        """

rule genome_align_no_plasmid:
    input:
        genome_bam=join(FULL_GENOME_ALIGN_DIR, '{sample}.{biorep}.bam'),
        plasmid_bam=join(FULL_PLASMID_ALIGN_DIR, '{sample}.{biorep}.bam')
    output:
        bam=join(GENOME_ALIGN_NOPLASMID_DIR, '{sample}.{biorep}.bam'),
        bai=join(GENOME_ALIGN_NOPLASMID_DIR, '{sample}.{biorep}.bam.bai')
    run:
        nanoporepipeline.remove_plasmid_alignments(input.genome_bam, input.plasmid_bam, output.bam)

rule extract_indels:
    input:
        bam=join(GENOME_ALIGN_NOPLASMID_DIR, '{sample}.{biorep}.bam')
    output:
        indels=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.indels.tsv')
    run:
        nanoporepipeline.extract_indels(input.bam, output.indels)

rule cluster_indels:
    input:
        indels=expand(join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.indels.tsv'), sample=samples_biorep)
    output:
        outtsv=join(RAW_STRUCTURAL_VARIANTS_DIR, 'clustered_indels.tsv')
    threads: 96
    run:
        nanoporepipeline.cluster_indels(input.indels, output.outtsv, threads)

rule filter_indels:
    input:
        tsv=join(RAW_STRUCTURAL_VARIANTS_DIR, 'clustered_indels.tsv')
    output:
        tsv1=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indels.tsv'),
        tsv2=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indels.reads.tsv')
    shell:
        """
        Rscript {NANOPOREPIPELINE_FILTER_INDELS_SCRIPT_PATH} {input.tsv} {output.tsv1} {output.tsv2}
        """

rule fgsv_svpileup:
    input:
        bam=join(GENOME_ALIGN_NOPLASMID_DIR, '{sample}.{biorep}.bam')
    output:
        outtxt=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.fgsv_svpileup.txt'),
        outbam=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.fgsv_svpileup.bam')
    params:
        prefix=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}')
    shell:
        """
        fgsv.jar SvPileup --slop 20 -i {input.bam} -o {params.prefix}
        mv {params.prefix}.txt {output.outtxt}
        mv {params.prefix}.bam {output.outbam}
        """

rule fgsv_filter_size:
    input:
        txt=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.fgsv_svpileup.txt'),
        genome=GENOME_FASTA
    output:
        outtxt=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.fgsv_svpileup.sizefilt.txt')
    run:
        nanoporepipeline.fgsv_filter_size(input.txt, input.genome, output.outtxt)

rule fgsv_aggsvpileup:
    input:
        svpileup=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.fgsv_svpileup.sizefilt.txt'),
        bam=join(FULL_GENOME_ALIGN_DIR, '{sample}.{biorep}.bam')
    output:
        outtsv=join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.{biorep}.aggsvpileup.tsv')
    shell:
        """
        fgsv.jar AggregateSvPileup -i {input.svpileup} -b {input.bam} --min-frequency 0.000001 --min-breakpoint-support 1 -d 20 -o {output.outtsv}
        """

rule merge_breakpoints:
    input:
        tsvs=expand(join(RAW_STRUCTURAL_VARIANTS_DIR, '{sample}.aggsvpileup.tsv'), sample=samples_biorep)
    output:
        tsv=join(RAW_STRUCTURAL_VARIANTS_DIR, 'merged_breakpoints.tsv')
    run:
        nanoporepipeline.merge_breakpoints(input.tsvs, output.tsv)

rule filter_breakpoints:
    input:
        tsv=join(RAW_STRUCTURAL_VARIANTS_DIR, 'merged_breakpoints.tsv')
    output:
        tsv1=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'breakpoints.tsv')
    shell:
        """
        Rscript {NANOPOREPIPELINE_FILTER_BREAKPOINTS_SCRIPT_PATH} {input.tsv} {output.tsv1}
        """

rule overlapping_breakpoints_indels:
    input:
        indels=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indels.tsv'),
        breakpoints=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'breakpoints.tsv')
    output:
        tsv=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'overlapping_breakpoints_indels.tsv')
    run:
        nanoporepipeline.overlapping_breakpoints_indels(input.indels, input.breakpoints, output.tsv)

rule overlapping_insertion_sites:
    input:
        indels=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indels.tsv'),
        breakpoints=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'breakpoints.tsv'),
        tsv1=join(INSERTION_ANALYSIS_DIR, 'programmed_insertion_sites.tsv'),
        tsv2=join(INSERTION_ANALYSIS_DIR, 'wt_insertion_sites.tsv')
    output:
        tsv=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'overlapping_insertion_sites.tsv')
    run:
        nanoporepipeline.overlapping_insertion_sites(input.indels, input.breakpoints, input.tsv1, input.tsv2, output.tsv)

rule indel_site_coverage:
    input:
        indels=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indels.tsv'),
        bam=expand(join(FULL_GENOME_ALIGN_DIR, '{sample}.bam'), sample=samples_biorep)
    output:
        tsv=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indel_site_coverage.tsv')
    threads: 16
    run:
        nanoporepipeline.indel_site_coverage(input.indels, input.bam, threads, output.tsv)

rule final_strucvars:
    input:
        indels=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indels.tsv'),
        indels_cov=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'indel_site_coverage.tsv'),
        breakpoints=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'breakpoints.tsv'),
        overlapping_breakpoints_indels=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'overlapping_breakpoints_indels.tsv'),
        overlapping_insertion_sites=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'overlapping_insertion_sites.tsv')
    output:
        indels=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'final_indels.tsv'),
        breakpoints=join(FINAL_STRUCTURAL_VARIANTS_DIR, 'final_breakpoints.tsv')
    shell:
        """
        Rscript {NANOPOREPIPELINE_FINAL_STRUCVARS_SCRIPT_PATH} {input.indels} {input.indels_cov} {input.breakpoints} \
          {input.overlapping_breakpoints_indels} {input.overlapping_insertion_sites} {output.indels} {output.breakpoints}
        """
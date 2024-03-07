from os.path import join, dirname
import sys, pkg_resources

BLASTP_PATH = 'blastp'

LINEARFOLD_PATH = pkg_resources.resource_filename('bridgerna2024', 'bin/linux/LinearFold-1.0/linearfold')
BBDUK_PATH = pkg_resources.resource_filename('bridgerna2024', 'bin/linux/bbmap/bbduk.sh')
BBDUK_ADAPTERS_PATH = pkg_resources.resource_filename('bridgerna2024', 'bin/linux/bbmap/resources/adapters.fa')
BWA_PATH = pkg_resources.resource_filename('bridgerna2024', 'bin/linux/bwa')

RNASEQ_SNAKEFILE_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/rnaseq.Snakefile')
RNASEQ_SNAKEMAKE_CONFIG_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/rnaseq.config.yml')
DONORSCREEN_SNAKEFILE_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/donorscreen.Snakefile')
DONORSCREEN_SNAKEMAKE_CONFIG_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/donorscreen.config.yml')
TARGETSCREEN_SNAKEFILE_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/targetscreen.Snakefile')
TARGETSCREEN_SNAKEMAKE_CONFIG_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/targetscreen.config.yml')
NANOPOREPIPELINE_SNAKEFILE_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/nanoporepipeline.Snakefile')
NANOPOREPIPELINE_SNAKEMAKE_CONFIG_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/nanoporepipeline.config.yml')
STRUCTURE_SNAKEFILE_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/structure.Snakefile')
STRUCTURE_SNAKEMAKE_CONFIG_PATH = pkg_resources.resource_filename('bridgerna2024', 'snakemake/structure.config.yml')

MMSEQS2_PATH = 'mmseqs'
MUSCLE_PATH = 'muscle'
MAFFT_QINSI_PATH = 'mafft-qinsi'
MAFFT_XINSI_PATH = 'mafft-xinsi'

# Script paths
RNASEQ_VISUALIZE_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/rnaseq_visualize.R')
STRUCTURE_VISUALIZE_ALIGNMENT_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/structure_visualize_alignment.R')
TARGETSCREEN_OLIGO_COVERAGE_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/targetscreen_oligo_coverage.R')
TARGETSCREEN_VISUALIZE_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/targetscreen_visualize.R')
DONORSCREEN_VISUALIZE_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/donorscreen_visualize.R')
NANOPOREPIPELINE_FILTER_BREAKPOINTS_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/nanoporepipeline_filter_breakpoints.R')
NANOPOREPIPELINE_FILTER_INDELS_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/nanoporepipeline_filter_indels.R')
NANOPOREPIPELINE_FINAL_STRUCVARS_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/nanoporepipeline_final_strucvars.R')
NANOPOREPIPELINE_INSERTION_ANALYSIS_SCRIPT_PATH = pkg_resources.resource_filename('bridgerna2024', 'scripts/nanoporepipeline_insertion_analysis.R')
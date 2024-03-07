import click
import os, random, sys
random.seed(42)

from Bio import SeqIO

from bridgerna2024 import structure
from bridgerna2024 import rnaseq
from bridgerna2024 import donorscreen
from bridgerna2024 import targetscreen
from bridgerna2024 import nanoporepipeline
from bridgerna2024.help import CustomHelp

@click.group(cls=CustomHelp)
def cli():
    """Code repository accompanying Durrant & Perry et al. 2024"""
    pass

@cli.command('rnaseq', help_priority=1)
@click.argument('workdir')
@click.option("--threads", "-t", default=1)
def rnaseq_command(workdir, threads):
    """Run RNA-seq snakemake pipeline."""
    rnaseq._run_rnaseq_snakemake(workdir, threads)

@cli.command('targetscreen', help_priority=3)
@click.argument('workdir')
@click.option("--threads", "-t", default=1)
def targetscreen_command(workdir, threads):
    """Run target screen snakemake pipeline."""
    targetscreen._run_targetscreen_snakemake(workdir, threads)

@cli.command('donorscreen', help_priority=4)
@click.argument('workdir')
@click.option("--threads", "-t", default=1)
def donorscreen_command(workdir, threads):
    """Run donor screen snakemake pipeline."""
    donorscreen._run_donorscreen_snakemake(workdir, threads)

@cli.command('nanoporepipeline', help_priority=5)
@click.argument('workdir')
@click.option("--threads", "-t", default=1)
def nanoporepipeline_command(workdir, threads):
    """Run genome insertion nanopore snakemake pipeline."""
    nanoporepipeline._run_nanoporepipeline_snakemake(workdir, threads)

@cli.group('structure', cls=CustomHelp)
def structure_group():
    """Run bridge RNA structure analyses."""
    pass

@structure_group.command('workflow', short_help="Run the full bridgeRNA structure discovery pipeline.", help_priority=1)
@click.argument("workdir")
@click.option("--threads", "-t", default=1)
@click.option("--is110-db-path", "-idp", default=None)
@click.option("--is110-blastdb-path", "-ibp", default=None)
@click.option("--retrieve-min-pident", default=None)
@click.option("--retrieve-min-alnlen", default=None)
@click.option("--retrieve-max-unique-proteins", default=None)
@click.option("--retrieve-max-flank-length", default=None)
@click.option("--retrieve-fiveprime-flank-offset", default=None)
@click.option("--retrieve-threeprime-flank-offset", default=None)
@click.option("--filter-fiveprime-flank-length", default=None)
@click.option("--filter-threeprime-flank-length", default=None)
@click.option("--filter-max-fiveprime-flank-length", default=None)
@click.option("--filter-max-threeprime-flank-length", default=None)
@click.option("--filter-cluster-tpase-min-pident", default=None)
@click.option("--filter-cluster-tpase-min-alnlen", default=None)
@click.option("--filter-cluster-flank-min-pident", default=None)
@click.option("--filter-cluster-flank-min-alnlen", default=None)
@click.option("--filter-infer-flank-length", default=None)
@click.option("--filter-max-seqs", default=None)
@click.option("--filter-random-seq", default=None)
@click.option("--foldalign-aligner", default=None)
@click.option("--foldalign-max-gaps", default=None)
@click.option("--foldalign-fiveprime-user-seq", default=None)
@click.option("--foldalign-threeprime-user-seq", default=None)
@click.option("--foldalign-revcomp-fold", default=None)
@click.option("--foldalign-keep-user-seq-cols", default=None)
def structure_workflow_command(
        workdir, threads, is110_db_path=None, is110_blastdb_path=None,
        retrieve_min_pident=None, retrieve_min_alnlen=None,
        retrieve_max_unique_proteins=None, retrieve_max_flank_length=None, retrieve_fiveprime_flank_offset=None,
        retrieve_threeprime_flank_offset=None, filter_fiveprime_flank_length=None, filter_threeprime_flank_length=None,
        filter_max_fiveprime_flank_length=None, filter_max_threeprime_flank_length=None,
        filter_cluster_tpase_min_pident=None, filter_cluster_tpase_min_alnlen=None,
        filter_cluster_flank_min_pident=None, filter_cluster_flank_min_alnlen=None, filter_infer_flank_length=None,
        filter_max_seqs=None, filter_random_seq=None, foldalign_aligner=None, foldalign_max_gaps=None,
        foldalign_fiveprime_user_seq=None, foldalign_threeprime_user_seq=None, foldalign_revcomp_fold=None,
        foldalign_keep_user_seq_cols=None
):

    """Run the full bridgeRNA structure prediction workflow search as a snakemake pipeline."""
    structure._run_structure_workflow(
        workdir, threads, is110_db_path, is110_blastdb_path, retrieve_min_pident, retrieve_min_alnlen,
        retrieve_max_unique_proteins, retrieve_max_flank_length, retrieve_fiveprime_flank_offset,
        retrieve_threeprime_flank_offset, filter_fiveprime_flank_length, filter_threeprime_flank_length,
        filter_max_fiveprime_flank_length, filter_max_threeprime_flank_length,
        filter_cluster_tpase_min_pident, filter_cluster_tpase_min_alnlen,
        filter_cluster_flank_min_pident, filter_cluster_flank_min_alnlen, filter_infer_flank_length,
        filter_max_seqs, filter_random_seq, foldalign_aligner, foldalign_max_gaps,
        foldalign_fiveprime_user_seq, foldalign_threeprime_user_seq, foldalign_revcomp_fold,
        foldalign_keep_user_seq_cols
    )

@structure_group.command('search', short_help="Search for related IS110s.", help_priority=2)
@click.argument("fasta")
@click.option("--is110-blastdb-path", "-ibp", default='WORKDIR/00.database/IS110_blastdb/IS110_tpases.faa')
@click.option("--outdir", "-o", default='structure_output')
@click.option("--prefix", "-p", default='01')
@click.option("--threads", "-t", default=1)
@click.option("--force/--no-force", default=False)
def structure_search_command(fasta, is110_blastdb_path, outdir, prefix, threads, force):
    """Search IS110 protein sequence in FASTA against all known IS110 proteins."""
    if os.path.isdir(outdir):
        if not force:
            print("Output directory {}/ already exists, use --force to overwrite. Exiting...".format(outdir))
            sys.exit()
        else:
            print("Using existing output directory {}/...".format(outdir))

    query_seq = None
    for rec in SeqIO.parse(fasta, 'fasta'):
        query_seq = {'name': rec.id, 'seq': str(rec.seq)}
        break

    structure._search(query_seq, is110_blastdb_path, outdir, prefix, threads)


@structure_group.command('retrieve', short_help="Filter protein sequences by % identity and retrieve CDS flanks.", help_priority=3)
@click.argument("blastp")
@click.option("--is110-db-path", "-idp", default='WORKDIR/00.database/IS110.db')
@click.option("--min-pident", "-mpi", default=30.0)
@click.option("--min-alnlen", "-mal", default=80.0)
@click.option("--max-unique-proteins", "-mup", default=2000)
@click.option("--max-flank-length", "-mfl", default=2000)
@click.option("--fiveprime-flank-offset", "-5fo", default=0)
@click.option("--threeprime-flank-offset", "-3fo", default=0)
@click.option("--outdir", "-o", default='structure_output')
@click.option("--prefix", "-p", default='02')
@click.option("--force/--no-force", default=False)
def structure_retrieve_command(blastp, is110_db_path, min_pident, min_alnlen, max_unique_proteins, max_flank_length, fiveprime_flank_offset, threeprime_flank_offset, outdir, prefix, force):
    if os.path.isdir(outdir):
        if not force:
            print("Output directory {}/ already exists, use --force to overwrite. Exiting...".format(outdir))
            sys.exit()
        else:
            print("Using existing output directory {}/...".format(outdir))

    structure._retrieve(blastp, is110_db_path, min_pident, min_alnlen, max_unique_proteins, max_flank_length,
                        fiveprime_flank_offset, threeprime_flank_offset, outdir, prefix)


@structure_group.command('filter', short_help="Filter sequences by clustering.", help_priority=4)
@click.argument("retrieved")
@click.option("--fiveprime-flank-length", "-5fl", default=300)
@click.option("--threeprime-flank-length", "-3fl", default=300)
@click.option("--max-fiveprime-flank-length", "-5fl", default=500)
@click.option("--max-threeprime-flank-length", "-3fl", default=500)
@click.option("--cluster-tpase-min-pident", "-ctid", default=98.0)
@click.option("--cluster-tpase-min-alnlen", "-ctal", default=80.0)
@click.option("--cluster-flank-min-pident", "-cfid", default=95.0)
@click.option("--cluster-flank-min-alnlen", "-cfal", default=80.0)
@click.option("--infer-flank-ength/--no-infer-flank-length", default=False)
@click.option("--max-seqs", "-ms", default=200)
@click.option("--random-seq/--closest-seq", default=True)
@click.option("--outdir", "-o", default='structure_output')
@click.option("--prefix", "-p", default='03')
@click.option("--threads", "-t", default=1)
@click.option("--force/--no-force", default=False)
def structure_filter_command(retrieved, fiveprime_flank_length, threeprime_flank_length, max_fiveprime_flank_length,
                             max_threeprime_flank_length, cluster_tpase_min_pident, cluster_tpase_min_alnlen,
                             cluster_flank_min_pident, cluster_flank_min_alnlen, infer_flank_length, max_seqs,
                             random_seq, outdir, prefix, threads, force):

    if os.path.isdir(outdir):
        if not force:
            print("Output directory {}/ already exists, use --force to overwrite. Exiting...".format(outdir))
            sys.exit()
        else:
            print("Using existing output directory {}/...".format(outdir))

    structure._filter(retrieved, fiveprime_flank_length, threeprime_flank_length, max_fiveprime_flank_length,
                      max_threeprime_flank_length, cluster_tpase_min_pident, cluster_tpase_min_alnlen,
                      cluster_flank_min_pident, cluster_flank_min_alnlen, infer_flank_length, max_seqs, random_seq,
                      outdir, prefix, threads)

@structure_group.command('fold-align', short_help="Fold dRNA sequences using linearfold and align with muscle.", help_priority=5)
@click.argument("filtered")
@click.option("--aligner", "-a", type=click.Choice(['muscle', 'mafft-xinsi', 'mafft-qinsi']), default='muscle')
@click.option("--max-gaps", "-mg", default=50.0)
@click.option("--fiveprime-user-seq", "-5us", default=None)
@click.option("--threeprime-user-seq", "-3us", default=None)
@click.option("--revcomp-fold/--no-revcomp-fold", default=False)
@click.option("--keep-user-seq-cols/--no-keep-user-seq-cols", default=False)
@click.option("--outdir", "-o", default='structure_output')
@click.option("--prefix", "-p", default='04')
@click.option("--threads", "-t", default=1)
@click.option("--force/--no-force", default=False)
def structure_fold_align_command(filtered, aligner, max_gaps, fiveprime_user_seq, threeprime_user_seq, revcomp_fold,
                                 keep_user_seq_cols, outdir, prefix, threads, force):

    if os.path.isdir(outdir):
        if not force:
            print("Output directory {}/ already exists, use --force to overwrite. Exiting...".format(outdir))
            sys.exit()
        else:
            print("Using existing output directory {}/...".format(outdir))

    if fiveprime_user_seq is not None and threeprime_user_seq is None:
        print("If specifying the user sequences, you must specify both the --fiveprime-user-seq and the --threeprime-user-seq")
        return
    elif fiveprime_user_seq is None and threeprime_user_seq is not None:
        print("If specifying the user sequences, you must specify both the --fiveprime-user-seq and the --threeprime-user-seq")
        return

    structure._fold_align(filtered, aligner, max_gaps, fiveprime_user_seq, threeprime_user_seq, outdir, prefix, threads, revcomp_fold, keep_user_seq_cols, force)
from glob import glob
from os.path import basename, join, isfile

from Bio import SeqIO

from bridgerna2024 import *
from bridgerna2024 import structure

from collections import defaultdict

WD = config['workdir']
DATABASE_DIR = join(WD, config['database_dir'])
QUERY_DIR = join(WD, config['query_dir'])
SEARCH_DIR = join(WD, config['search_dir'])
RETRIEVE_DIR = join(WD, config['retrieve_dir'])
FILTER_DIR = join(WD, config['filter_dir'])
FOLDALIGN_DIR = join(WD, config['foldalign_dir'])
VISUALIZE_DIR = join(WD, config['visualize_dir'])

IS110_DB = join(DATABASE_DIR, config['is110_db'])
IS110_TPASE_BLASTDB = join(DATABASE_DIR, config['is110_blastdb'])
QUERIES = list(set([basename(f).split('.')[0] for f in glob(join(QUERY_DIR, '*fasta'))]))

QUERY_METADATA = join(WD, "metadata.tsv")
QUERY_FIVEPRIME_FLANK_LENGTH = dict()
QUERY_THREEPRIME_FLANK_LENGTH = dict()

if isfile(QUERY_METADATA):

    with open(QUERY_METADATA, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip().split('\t')
            QUERY_FIVEPRIME_FLANK_LENGTH[line[0]] = int(line[1])
            QUERY_THREEPRIME_FLANK_LENGTH[line[0]] = int(line[2])

rule all:
    input:
        [join(VISUALIZE_DIR, '{query}.png'.format(query=qry)) for qry in QUERIES]

rule search:
    input:
        query = join(QUERY_DIR, '{query}.fasta'),
    output:
        blastp = join(SEARCH_DIR, '{query}.search.blastp.tsv')
    params:
        outdir = SEARCH_DIR,
        prefix = '{query}'
    threads: 8
    run:
        query_seq = None
        for rec in SeqIO.parse(input.query, 'fasta'):
            query_seq = {'name': rec.id, 'seq': str(rec.seq)}
            break
        structure._search(
            query_seq, IS110_TPASE_BLASTDB, params.outdir, params.prefix, threads
        )

rule retrieve:
    input:
        blastp = join(SEARCH_DIR, '{query}.search.blastp.tsv')
    output:
        retrieved = join(RETRIEVE_DIR, '{query}.retrieve.tsv')
    params:
        outdir = RETRIEVE_DIR,
        prefix = '{query}',
        min_pident = config['retrieve_min_pident'],
        min_alnlen = config['retrieve_min_alnlen'],
        max_unique_proteins = config['retrieve_max_unique_proteins'],
        max_flank_length = config['retrieve_max_flank_length'],
        fiveprime_flank_offset = config['retrieve_fiveprime_flank_offset'],
        threeprime_flank_offset = config['retrieve_threeprime_flank_offset']
    run:
        structure._retrieve(
            input.blastp, IS110_DB, params.min_pident, params.min_alnlen, params.max_unique_proteins,
            params.max_flank_length, params.fiveprime_flank_offset, params.threeprime_flank_offset, params.outdir,
            params.prefix
        )

rule filter:
    input:
        retrieved = join(RETRIEVE_DIR, '{query}.retrieve.tsv')
    output:
        filtered = join(FILTER_DIR, '{query}.filtered.tsv')
    params:
        outdir = FILTER_DIR,
        prefix = '{query}',
        fiveprime_flank_length = config['filter_fiveprime_flank_length'],
        threeprime_flank_length = config['filter_threeprime_flank_length'],
        max_fiveprime_flank_length = config['filter_max_fiveprime_flank_length'],
        max_threeprime_flank_length = config['filter_max_threeprime_flank_length'],
        cluster_tpase_min_pident = config['filter_cluster_tpase_min_pident'],
        cluster_tpase_min_alnlen = config['filter_cluster_tpase_min_alnlen'],
        cluster_flank_min_pident = config['filter_cluster_flank_min_pident'],
        cluster_flank_min_alnlen = config['filter_cluster_flank_min_alnlen'],
        infer_flank_length = config['filter_infer_flank_length'],
        max_seqs = config['filter_max_seqs'],
        random_seq = config['filter_random_seq']
    threads: 8
    run:
        if len(QUERY_FIVEPRIME_FLANK_LENGTH) == 0:
            structure._filter(
                input.retrieved, params.fiveprime_flank_length, params.threeprime_flank_length,
                params.max_fiveprime_flank_length, params.max_threeprime_flank_length,
                params.cluster_tpase_min_pident, params.cluster_tpase_min_alnlen, params.cluster_flank_min_pident,
                params.cluster_flank_min_alnlen, params.infer_flank_length, params.max_seqs, params.random_seq,
                params.outdir, params.prefix, threads
            )
        else:
            structure._filter(
                input.retrieved, QUERY_FIVEPRIME_FLANK_LENGTH[params.prefix], QUERY_THREEPRIME_FLANK_LENGTH[params.prefix],
                params.max_fiveprime_flank_length, params.max_threeprime_flank_length,
                params.cluster_tpase_min_pident, params.cluster_tpase_min_alnlen, params.cluster_flank_min_pident,
                params.cluster_flank_min_alnlen, False, params.max_seqs, params.random_seq,
                params.outdir, params.prefix, threads
            )

rule foldalign:
    input:
        filtered = join(FILTER_DIR, '{query}.filtered.tsv')
    output:
        foldalign = join(FOLDALIGN_DIR, '{query}.fold_align.tsv')
    params:
        outdir = FOLDALIGN_DIR,
        prefix = '{query}',
        aligner = config['foldalign_aligner'],
        max_gaps = config['foldalign_max_gaps'],
        fiveprime_user_seq = config['foldalign_fiveprime_user_seq'],
        threeprime_user_seq = config['foldalign_threeprime_user_seq'],
        revcomp_fold = config['foldalign_revcomp_fold'],
        keep_user_seq_cols = config['foldalign_keep_user_seq_cols']
    threads: 8
    run:
        structure._fold_align(
            input.filtered, params.aligner, params.max_gaps, params.fiveprime_user_seq,
            params.threeprime_user_seq, params.outdir, params.prefix, threads,
            params.revcomp_fold, params.keep_user_seq_cols, force=False
        )

rule visualize:
    input:
        foldalign = join(FOLDALIGN_DIR, '{query}.fold_align.tsv')
    output:
        viz = join(VISUALIZE_DIR, '{query}.png')
    params:
        outdir = VISUALIZE_DIR,
        prefix = '{query}',
        fiveprime_flank_offset = config['retrieve_fiveprime_flank_offset'],
        threeprime_flank_offset = config['retrieve_threeprime_flank_offset']
    shell:
        """
        Rscript {STRUCTURE_VISUALIZE_ALIGNMENT_SCRIPT_PATH} {input.foldalign} \
        {params.fiveprime_flank_offset} {params.threeprime_flank_offset} \
        {params.prefix} {params.outdir}
        """
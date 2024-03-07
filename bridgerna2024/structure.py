import os
from Bio import SeqIO
from Bio.Seq import reverse_complement
import shutil
from bridgerna2024 import *
from bridgerna2024 import misc

import subprocess
import pandas as pd
from collections import defaultdict
import random
import sqlite3
from tqdm import tqdm
from multiprocessing import Pool
import numpy as np
from forgi.graph import bulge_graph


def _run_structure_workflow(
        workdir, threads, is110_db_path=None, is110_blastdb_path=None, retrieve_min_pident=None, retrieve_min_alnlen=None,
        retrieve_max_unique_proteins=None, retrieve_max_flank_length=None, retrieve_fiveprime_flank_offset=None,
        retrieve_threeprime_flank_offset=None, filter_fiveprime_flank_length=None, filter_threeprime_flank_length=None,
        filter_max_fiveprime_flank_length=None, filter_max_threeprime_flank_length=None,
        filter_cluster_tpase_min_pident=None, filter_cluster_tpase_min_alnlen=None,
        filter_cluster_flank_min_pident=None, filter_cluster_flank_min_alnlen=None, filter_infer_flank_length=None,
        filter_max_seqs=None, filter_random_seq=None, foldalign_aligner=None, foldalign_max_gaps=None,
        foldalign_fiveprime_user_seq=None, foldalign_threeprime_user_seq=None, foldalign_revcomp_fold=None,
        foldalign_keep_user_seq_cols=None
):

    cmd = f"snakemake --snakefile {STRUCTURE_SNAKEFILE_PATH} --configfile {STRUCTURE_SNAKEMAKE_CONFIG_PATH} -j {threads} " \
          f"--config workdir={workdir} "

    if is110_db_path is not None:
        cmd += f"is110_db_path={is110_db_path} "
    if is110_blastdb_path is not None:
        cmd += f"is110_blastdb_path={is110_blastdb_path} "
    if retrieve_min_pident is not None:
        cmd += f"retrieve_min_pident={retrieve_min_pident} "
    if retrieve_min_alnlen is not None:
        cmd += f"retrieve_min_alnlen={retrieve_min_alnlen} "
    if retrieve_max_unique_proteins is not None:
        cmd += f"retrieve_max_unique_proteins={retrieve_max_unique_proteins} "
    if retrieve_max_flank_length is not None:
        cmd += f"retrieve_max_flank_length={retrieve_max_flank_length} "
    if retrieve_fiveprime_flank_offset is not None:
        cmd += f"retrieve_fiveprime_flank_offset={retrieve_fiveprime_flank_offset} "
    if retrieve_threeprime_flank_offset is not None:
        cmd += f"retrieve_threeprime_flank_offset={retrieve_threeprime_flank_offset} "
    if filter_fiveprime_flank_length is not None:
        cmd += f"filter_fiveprime_flank_length={filter_fiveprime_flank_length} "
    if filter_threeprime_flank_length is not None:
        cmd += f"filter_threeprime_flank_length={filter_threeprime_flank_length} "
    if filter_max_fiveprime_flank_length is not None:
        cmd += f"filter_max_fiveprime_flank_length={filter_max_fiveprime_flank_length} "
    if filter_max_threeprime_flank_length is not None:
        cmd += f"filter_max_threeprime_flank_length={filter_max_threeprime_flank_length} "
    if filter_cluster_tpase_min_pident is not None:
        cmd += f"filter_cluster_tpase_min_pident={filter_cluster_tpase_min_pident} "
    if filter_cluster_tpase_min_alnlen is not None:
        cmd += f"filter_cluster_tpase_min_alnlen={filter_cluster_tpase_min_alnlen} "
    if filter_cluster_flank_min_pident is not None:
        cmd += f"filter_cluster_flank_min_pident={filter_cluster_flank_min_pident} "
    if filter_cluster_flank_min_alnlen is not None:
        cmd += f"filter_cluster_flank_min_alnlen={filter_cluster_flank_min_alnlen} "
    if filter_infer_flank_length is not None:
        cmd += f"filter_infer_flank_length={filter_infer_flank_length} "
    if filter_max_seqs is not None:
        cmd += f"filter_max_seqs={filter_max_seqs} "
    if filter_random_seq is not None:
        cmd += f"filter_random_seq={filter_random_seq} "
    if foldalign_aligner is not None:
        cmd += f"foldalign_aligner={foldalign_aligner} "
    if foldalign_max_gaps is not None:
        cmd += f"foldalign_max_gaps={foldalign_max_gaps} "
    if foldalign_fiveprime_user_seq is not None:
        cmd += f"foldalign_fiveprime_user_seq={foldalign_fiveprime_user_seq} "
    if foldalign_threeprime_user_seq is not None:
        cmd += f"foldalign_threeprime_user_seq={foldalign_threeprime_user_seq} "
    if foldalign_revcomp_fold is not None:
        cmd += f"foldalign_revcomp_fold={foldalign_revcomp_fold} "
    if foldalign_keep_user_seq_cols is not None:
        cmd += f"foldalign_keep_user_seq_cols={foldalign_keep_user_seq_cols} "

    """Run the full bridgeRNA structure prediction workflow search as a snakemake pipeline."""
    print("Running the structure prediction workflow...")
    print(cmd)
    subprocess.run(cmd.split())

def _search(query_seq, is110_blastdb_path, outdir, prefix, threads):

    print("Running search...")
    misc.print_variables_dictionary(locals())

    os.makedirs(outdir, exist_ok=True)
    tmpdir = os.path.join(outdir, 'tmp_' + str(random.randint(-sys.maxsize, sys.maxsize)))
    os.makedirs(tmpdir, exist_ok=True)

    query_faa = os.path.join(tmpdir, 'query.faa')
    with open(query_faa, 'w') as outfile:
        print('>' + query_seq['name'], query_seq['seq'], sep='\n', file=outfile)

    # Running blastp of protein sequences
    cmd = '{blastp} -query {query} -db {db} -num_threads {t} -max_target_seqs 1000000 -evalue 1e-6 -outfmt'.format(
        blastp=BLASTP_PATH, query=query_faa, db=is110_blastdb_path, t=threads
    )
    outfmt = "6 qseqid sseqid qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    blast_out_path = os.path.join(tmpdir, 'query.proteins.blastp.tsv')
    outfile = open(blast_out_path, 'w')
    print("Executing command:", cmd + " " + outfmt)
    subprocess.run(cmd.split() + [outfmt], stdout=outfile)
    outfile.close()

    # Parsing results
    blastp_results = pd.read_csv(blast_out_path, sep='\t', header=None)
    blastp_results.columns = ['qseqid', 'sseqid', 'qlen', 'slen', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                              'qend', 'sstart', 'send', 'evalue', 'bitscore']

    blastp_out = open(os.path.join(outdir, prefix + '.search.blastp.tsv'), 'w')
    print('p100_id', 'alnlen', 'qlen', 'plen', 'pident', 'evalue', 'bitscore', sep='\t', file=blastp_out)

    print("Processing results...")
    keep_proteins = set()
    for sseqid, length, qlen, slen, pident, evalue, bitscore in zip(blastp_results.sseqid, blastp_results.length,
                                                                    blastp_results.qlen, blastp_results.slen,
                                                                    blastp_results.pident, blastp_results.evalue,
                                                                    blastp_results.bitscore):
        print(sseqid, length, qlen, slen, pident, evalue, bitscore, sep='\t', file=blastp_out)

    blastp_out.close()

    shutil.copyfile(query_faa, os.path.join(outdir, prefix + '.search.query.faa'))

    shutil.rmtree(tmpdir)

    print("Search complete.")


def _retrieve(blastp_tsv, is110_db_path, min_pident, min_alnlen, max_unique_proteins, max_flank_length, fiveprime_flank_offset,
              threeprime_flank_offset, outdir, prefix):

    print("Running retrieve...")
    misc.print_variables_dictionary(locals())

    os.makedirs(outdir, exist_ok=True)
    tmpdir = os.path.join(outdir, 'tmp_' + str(random.randint(-sys.maxsize, sys.maxsize)))
    os.makedirs(tmpdir, exist_ok=True)

    keep_loci = set()

    blastp = pd.read_csv(blastp_tsv, sep='\t')
    blastp = blastp[blastp.pident > min_pident]
    blastp = blastp[(blastp.alnlen / blastp.qlen) * 100 > min_alnlen]
    blastp = blastp[(blastp.alnlen / blastp.plen) * 100 > min_alnlen]
    blastp = blastp.sort_values('pident', ascending=False)

    keep_proteins = set()
    for prot in blastp.p100_id:
        if len(keep_proteins) >= max_unique_proteins:
            break
        keep_proteins.add(prot)

    if len(keep_proteins) == 0:
        print("No proteins met filtering thresholds.")
        return False
    else:
        print("{} proteins met filtering thresholds...".format(len(keep_proteins)))

    p100_to_pident = dict()
    already = set()
    for p100, pident in zip(blastp.p100_id, blastp.pident):
        if p100 in already:
            continue
        already.add(p100)
        p100_to_pident[p100] = pident

    # Get protein flanks
    protein_flanks = []
    print(is110_db_path)
    with sqlite3.connect(is110_db_path) as conn:
        cur = conn.cursor()

        for prot in tqdm(keep_proteins, desc='Protein flanks retrieved'):
            cmd = f'SELECT seq FROM protein WHERE p100_id="{prot}"'
            prot_seq = None
            for res in cur.execute(cmd):
                prot_seq = res[0]

            cmd = f'SELECT contig_id,cds_start,cds_end,strand FROM coords WHERE p100_id="{prot}"'
            loci = []
            for res in cur.execute(cmd):
                loci.append(res)
            for contig_id, cds_start, cds_end, strand in loci:
                cmd = f'SELECT seq FROM contig WHERE contig_id="{contig_id}"'
                contig_seq = None
                for res in cur.execute(cmd):
                    contig_seq = res[0]
                    break

                if strand == '+':
                    cds = contig_seq[cds_start:cds_end]
                    cds_start += fiveprime_flank_offset
                    cds_end += threeprime_flank_offset
                else:
                    cds = reverse_complement(contig_seq[cds_start:cds_end])
                    cds_start -= threeprime_flank_offset
                    cds_end -= fiveprime_flank_offset

                leftflank = None
                if strand == '+':
                    start = cds_start - max_flank_length
                    end = cds_start
                    if start < 0:
                        start = 0
                    if end > len(contig_seq):
                        end = len(contig_seq)
                    leftflank = contig_seq[start:end]
                else:
                    start = cds_end
                    end = cds_end + max_flank_length
                    if start < 0:
                        start = 0
                    if end > len(contig_seq):
                        end = len(contig_seq)
                    leftflank = reverse_complement(contig_seq[start:end])

                rightflank = None
                if strand == '+':
                    start = cds_end
                    end = cds_end + max_flank_length
                    if start < 0:
                        start = 0
                    if end > len(contig_seq):
                        end = len(contig_seq)
                    rightflank = contig_seq[start:end]
                else:
                    start = cds_start - max_flank_length
                    end = cds_start
                    if start < 0:
                        start = 0
                    if end > len(contig_seq):
                        end = len(contig_seq)
                    rightflank = reverse_complement(contig_seq[start:end])

                protein_flanks.append(
                    (prot, contig_id, cds_start, cds_end, strand, leftflank, rightflank, cds, prot_seq))

        cur.close()

    try:
        protein_flanks = pd.DataFrame(protein_flanks)
        protein_flanks.columns = ['p100_id', 'contig_id', 'cds_start', 'cds_end', 'cds_strand', 'leftflank', 'rightflank',
                                      'cds', 'prot_seq']
    except ValueError:
        protein_flanks = pd.DataFrame(columns=['p100_id', 'contig_id', 'cds_start', 'cds_end', 'cds_strand', 'leftflank',
                                                  'rightflank', 'cds', 'prot_seq'])

    protein_flanks['max_flank_length'] = max_flank_length
    protein_flanks['leftflank_length'] = protein_flanks.leftflank.str.len()
    protein_flanks['rightflank_length'] = protein_flanks.rightflank.str.len()
    protein_flanks['cds_length'] = protein_flanks.cds_end - protein_flanks.cds_start
    protein_flanks['pident'] = [p100_to_pident[p100] for p100 in protein_flanks.p100_id]

    protein_flanks = protein_flanks[
        ['p100_id', 'pident', 'contig_id', 'cds_start', 'cds_end', 'cds_strand', 'max_flank_length', 'leftflank_length',
         'rightflank_length', 'cds_length', 'leftflank', 'rightflank', 'cds', 'prot_seq']]

    protein_flanks.sort_values(['pident', 'p100_id', 'contig_id', 'cds_start', 'cds_end'],
                               ascending=[False, True, True, True, True]).to_csv(
        os.path.join(outdir, prefix + '.retrieve.tsv'), sep='\t', index=None)

    shutil.rmtree(tmpdir)


def run_linearfold(seq):
    cmd = f'echo "{seq}" | {LINEARFOLD_PATH}'
    linearfold_out = subprocess.run(cmd, shell=True, capture_output=True).stdout.decode('utf-8').split()
    linearfold_out[-1] = float(linearfold_out[-1].replace("(", "").replace(")", ""))
    return linearfold_out


def map_linfold_to_aln(aln, linfold):
    linfold_aln = ''
    linfold_pos = 0
    for i, c in enumerate(aln):
        if c == '-':
            linfold_aln += '-'
        else:
            linfold_aln += linfold[linfold_pos]
            linfold_pos += 1

    return linfold_aln


def remove_gaps(seq, gaps):
    outseq = ''
    for i, c in enumerate(seq):
        if i in gaps:
            continue
        outseq += c
    return outseq


def average_length_protein(proteins, protein_seqs, expected_length=325):
    if len(proteins) == 1:
        return list(proteins)[0]

    if len(proteins) == 2:
        closest = None
        closest_dist = sys.maxsize
        for prot in proteins:
            dist = abs(len(protein_seqs[prot]) - expected_length)
            if dist < closest_dist:
                closest = prot
                closest_dist = dist
        return closest

    prot_lengths = dict()
    for prot in proteins:
        prot_lengths[prot] = len(protein_seqs[prot])

    median_length = np.median(list(prot_lengths.values()))
    closest = None
    closest_dist = sys.maxsize
    for prot in proteins:
        dist = abs(len(protein_seqs[prot]) - median_length)
        if dist < closest_dist:
            closest = prot
            closest_dist = dist

    return closest


def infer_left_flank_length(retrieved_df, tmpdir, fiveprime_flank_length, max_fiveprime_flank_length, flank_padding=50):
    leftflanks = []
    for lf in retrieved_df.leftflank:
        if len(leftflanks) == 10:
            break
        lf_leader = lf[-1000:-970]
        if len(lf_leader) != 30:
            continue

        mindiff = 30
        for other_leader, other_lf in leftflanks:
            ed = misc.edit_distance(lf_leader, other_leader)
            if ed < mindiff:
                mindiff = ed
        if mindiff <= 10:
            continue
        leftflanks.append((lf_leader, lf[-1000:]))

    tmp_leftflanks_fasta = os.path.join(tmpdir, 'inferlength_leftflanks.fasta')
    with open(tmp_leftflanks_fasta, 'w') as f:
        for i, (lf_leader, lf) in enumerate(leftflanks):
            f.write(f'>leftflank{i}\n{lf}\n')

    if len(leftflanks) != 10:
        print(f"Couldn't infer length of left flank, using default value of fiveprime_flank_length={fiveprime_flank_length}")
        return fiveprime_flank_length

    tmp_leftflanks_aln_fasta = tmp_leftflanks_fasta.replace('.fasta', '.aln.fasta')
    cmd = f'{MUSCLE_PATH} -align {tmp_leftflanks_fasta} -output {tmp_leftflanks_aln_fasta}'
    print("Running cmd:", cmd)
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    alnseqs = []
    alnlength = None
    nuc_counts = defaultdict(lambda: defaultdict(int))
    for rec in SeqIO.parse(tmp_leftflanks_aln_fasta, 'fasta'):
        recseq = str(rec.seq)
        alnseqs.append(recseq)
        alnlength = len(recseq)
        for i, c in enumerate(recseq):
            nuc_counts[i][c] += 1

    all_gaps_pos = list()
    top_nuc_props = dict()
    for i in nuc_counts:
        top_nuc = max(nuc_counts[i], key=lambda x: nuc_counts[i][x])
        try:
            second_top_nuc = sorted(nuc_counts[i], key=lambda x: nuc_counts[i][x], reverse=True)[1]
        except IndexError:
            second_top_nuc = None

        if top_nuc == '-':
            top_nuc = second_top_nuc

        if nuc_counts[i]['-'] == 9:
            all_gaps_pos.append(i)

        top_nuc_prop = nuc_counts[i][top_nuc] / sum(nuc_counts[i].values())
        top_nuc_props[i] = top_nuc_prop

    trim_left = 0
    if all_gaps_pos[0] == 0:
        for i in range(1, len(all_gaps_pos)):
            if all_gaps_pos[i] - all_gaps_pos[i-1] != 1:
                trim_left = i
                break

    trim_right = 0
    if all_gaps_pos[-1] == alnlength - 1:
        for i in range(len(all_gaps_pos)-2, -1, -1):
            if all_gaps_pos[i+1] - all_gaps_pos[i] != 1:
                trim_right = len(all_gaps_pos) - i - 1
                break

    top_nuc_idx_list = []
    top_nuc_props_list = []
    for i in range(trim_left, alnlength - trim_right):
        top_nuc_idx_list.append(i)
        top_nuc_props_list.append(top_nuc_props[i])

    max_split_idx = 0
    max_split_diff = -sys.maxsize
    for i in range(1, len(top_nuc_props_list)):
        diff = np.mean(top_nuc_props_list[i:]) - np.mean(top_nuc_props_list[:i])
        if diff > max_split_diff:
            max_split_idx = i
            max_split_diff = diff

    predicted_flanklengths = []
    for seq in alnseqs:
        seqsplit1 = seq[:max_split_idx]
        seqsplit2 = seq[max_split_idx:]

        seqsplit1_trim = seqsplit1[trim_left:]
        seqsplit2_trim = seqsplit2
        if trim_right > 0:
            seqsplit2_trim = seqsplit2[:-trim_right]

        predicted_flanklengths.append(len(seqsplit2_trim.replace('-', '')))

    inferred_flanklength = int(np.median(predicted_flanklengths) + flank_padding)

    if inferred_flanklength > max_fiveprime_flank_length:
        return max_fiveprime_flank_length
    else:
        return inferred_flanklength

def infer_right_flank_length(retrieved_df, tmpdir, threeprime_flank_length, max_threeprime_flank_length, flank_padding=50):
    rightflanks = []
    for rf in retrieved_df.rightflank:
        if len(rightflanks) == 10:
            break
        rf_leader = rf[970:1000]
        if len(rf_leader) != 30:
            continue

        mindiff = 30
        for other_leader, other_rf in rightflanks:
            ed = misc.edit_distance(rf_leader, other_leader)
            if ed < mindiff:
                mindiff = ed
        if mindiff <= 10:
            continue
        rightflanks.append((rf_leader, rf[:1000]))

    tmp_rightflanks_fasta = os.path.join(tmpdir, 'inferlength_rightflanks.fasta')
    with open(tmp_rightflanks_fasta, 'w') as f:
        for i, (lf_leader, lf) in enumerate(rightflanks):
            f.write(f'>rightflank{i}\n{lf}\n')

    if len(rightflanks) != 10:
        print(f"Couldn't infer length of right flank, using default value of threeprime_flank_length={threeprime_flank_length}")
        return threeprime_flank_length

    tmp_rightflanks_aln_fasta = tmp_rightflanks_fasta.replace('.fasta', '.aln.fasta')
    cmd = f'{MUSCLE_PATH} -align {tmp_rightflanks_fasta} -output {tmp_rightflanks_aln_fasta}'
    print("Running cmd:", cmd)
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    alnseqs = []
    alnlength = None
    nuc_counts = defaultdict(lambda: defaultdict(int))
    for rec in SeqIO.parse(tmp_rightflanks_aln_fasta, 'fasta'):
        recseq = str(rec.seq)
        alnseqs.append(recseq)
        alnlength = len(recseq)
        for i, c in enumerate(recseq):
            nuc_counts[i][c] += 1

    all_gaps_pos = list()
    top_nuc_props = dict()
    for i in nuc_counts:
        top_nuc = max(nuc_counts[i], key=lambda x: nuc_counts[i][x])
        try:
            second_top_nuc = sorted(nuc_counts[i], key=lambda x: nuc_counts[i][x], reverse=True)[1]
        except IndexError:
            second_top_nuc = None

        if top_nuc == '-':
            top_nuc = second_top_nuc

        if nuc_counts[i]['-'] == 9:
            all_gaps_pos.append(i)

        top_nuc_prop = nuc_counts[i][top_nuc] / sum(nuc_counts[i].values())
        top_nuc_props[i] = top_nuc_prop

    trim_left = 0
    if all_gaps_pos[0] == 0:
        for i in range(1, len(all_gaps_pos)):
            if all_gaps_pos[i] - all_gaps_pos[i-1] != 1:
                trim_left = i
                break

    trim_right = 0
    if all_gaps_pos[-1] == alnlength - 1:
        for i in range(len(all_gaps_pos)-2, -1, -1):
            if all_gaps_pos[i+1] - all_gaps_pos[i] != 1:
                trim_right = len(all_gaps_pos) - i - 1
                break

    top_nuc_idx_list = []
    top_nuc_props_list = []
    for i in range(trim_left, alnlength - trim_right):
        top_nuc_idx_list.append(i)
        top_nuc_props_list.append(top_nuc_props[i])

    max_split_idx = 0
    max_split_diff = -sys.maxsize
    for i in range(1, len(top_nuc_props_list)):
        diff = np.mean(top_nuc_props_list[:i]) - np.mean(top_nuc_props_list[i:])
        if diff > max_split_diff:
            max_split_idx = i
            max_split_diff = diff

    predicted_flanklengths = []
    for seq in alnseqs:
        seqsplit1 = seq[:max_split_idx]
        seqsplit2 = seq[max_split_idx:]

        seqsplit1_trim = seqsplit1[trim_left:]

        predicted_flanklengths.append(len(seqsplit1_trim.replace('-', '')))

    inferred_flanklength = int(np.median(predicted_flanklengths) + flank_padding)

    if inferred_flanklength > max_threeprime_flank_length:
        return max_threeprime_flank_length
    else:
        return inferred_flanklength


def _filter(retrieved_tsv, fiveprime_flank_length, threeprime_flank_length, max_fiveprime_flank_length,
            max_threeprime_flank_length, cluster_tpase_min_pident, cluster_tpase_min_alnlen,
            cluster_flank_min_pident, cluster_flank_min_alnlen, infer_flank_length, max_seqs, random_seq,
            outdir, prefix, threads):

    print("Running filter...")
    misc.print_variables_dictionary(locals())

    os.makedirs(outdir, exist_ok=True)
    tmpdir = os.path.join(outdir, 'tmp_' + str(random.randint(-sys.maxsize, sys.maxsize)))
    os.makedirs(tmpdir, exist_ok=True)

    retrieved_df = pd.read_csv(retrieved_tsv, sep='\t')
    retrieved_df = retrieved_df[~retrieved_df.prot_seq.str.contains("X")].reset_index()

    if retrieved_df is None or len(retrieved_df) == 0 or retrieved_df.shape[0] == 0:
        print("No sequences retrieved, exiting...")
        outdf = pd.DataFrame(columns=['recid', 'tpase_pident', 'left_flankseq', 'right_flankseq'])
        outdf.to_csv(os.path.join(outdir, prefix + '.filtered.tsv'), sep='\t', index=False)
        return

    if infer_flank_length:
        print("Inferring length of left and right flanks by alignment...")
        fiveprime_flank_length = infer_left_flank_length(retrieved_df, tmpdir, fiveprime_flank_length, max_fiveprime_flank_length)
        threeprime_flank_length = infer_right_flank_length(retrieved_df, tmpdir, threeprime_flank_length, max_threeprime_flank_length)

    print("Flank lengths:")
    misc.print_variables_dictionary(
        {'fiveprime_flank_length':fiveprime_flank_length, 'threeprime_flank_length':threeprime_flank_length}
    )

    # Cluster tpases with mmseqs2
    print("Clustering and filtering by Tpases...")
    p100_to_protseq = dict()
    for p100, protseq in zip(retrieved_df.p100_id, retrieved_df.prot_seq):
        p100_to_protseq[p100] = protseq

    prot_out_path = os.path.join(tmpdir, 'tpases.faa')
    outfile = open(prot_out_path, "w")
    for p100 in p100_to_protseq:
        print(">" + p100, p100_to_protseq[p100], sep='\n', file=outfile)
    outfile.close()

    minseqs_cmd = "{mmseqs} easy-linclust --threads {threads} --cov-mode 0 --min-seq-id {min_seq_id} -c {min_aln_len} {fna} {prefix} {tmpdir}".format(
        mmseqs=MMSEQS2_PATH, threads=threads, min_seq_id=cluster_tpase_min_pident / 100.0,
        min_aln_len=cluster_tpase_min_alnlen / 100.0, fna=prot_out_path,
        prefix=prot_out_path.replace('.faa', '.mmseqs'), tmpdir=tmpdir
    )
    print(minseqs_cmd)
    subprocess.run(minseqs_cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cluster2prots = defaultdict(set)
    with open(prot_out_path.replace('.faa', '.mmseqs') + "_cluster.tsv") as infile:
        for line in infile:
            line = line.strip().split()
            cluster2prots[line[0]].add(line[1])

    keep_proteins = set()
    for clust in cluster2prots:
        keep_proteins.add(average_length_protein(cluster2prots[clust], p100_to_protseq))

    keep_indices = set()
    for i, p100 in enumerate(retrieved_df.p100_id):
        if p100 in keep_proteins:
            keep_indices.add(i)
            keep_proteins.remove(p100)

    retrieved_df = retrieved_df.iloc[list(keep_indices)]

    # Extracting subflanks
    flanks_to_recid = dict()
    recid_to_flanks = dict()
    flanks_to_pident = dict()

    for leftflank, rightflank, pident in zip(retrieved_df.leftflank, retrieved_df.rightflank, retrieved_df.pident):

        try:
            if len(leftflank) < fiveprime_flank_length - 35:
                continue
            if len(leftflank) < 35:
                continue
            if len(rightflank) < threeprime_flank_length - 35:
                continue
            if len(rightflank) < 35:
                continue
        except:
            continue

        if len(leftflank) < fiveprime_flank_length:
            left_subflank = leftflank
        else:
            left_subflank = leftflank[len(leftflank) - fiveprime_flank_length:]

        if left_subflank.count('A') + left_subflank.count('C') + left_subflank.count('G') + left_subflank.count(
                'T') != len(left_subflank):
            continue

        if len(rightflank) < threeprime_flank_length:
            right_subflank = rightflank
        else:
            right_subflank = rightflank[:threeprime_flank_length]

        if right_subflank.count('A') + right_subflank.count('C') + right_subflank.count('G') + right_subflank.count(
                'T') != len(right_subflank):
            continue

        flanks = (left_subflank, right_subflank)
        if flanks not in flanks_to_recid:
            flanks_to_recid[flanks] = len(flanks_to_recid)
            recid_to_flanks[flanks_to_recid[flanks]] = flanks
            flanks_to_pident[flanks] = pident

    flanks_out_path = os.path.join(tmpdir, 'flanks.fna')
    with open(flanks_out_path, 'w') as out:
        for flanks in flanks_to_recid:
            recid = str(flanks_to_recid[flanks])
            print(">" + recid, flanks[0] + flanks[1], sep='\n', file=out)

    # Cluster with mmseqs2
    print("Clustering flanks...")
    minseqs_cmd = "{mmseqs} easy-linclust --threads {threads} --cov-mode 0 --min-seq-id {min_seq_id} -c {min_aln_len} {fna} {prefix} {tmpdir}"
    cmd1 = minseqs_cmd.format(
        mmseqs=MMSEQS2_PATH, threads=threads, min_seq_id=cluster_flank_min_pident / 100.0,
        min_aln_len=cluster_flank_min_alnlen / 100.0,
        fna=flanks_out_path, prefix=flanks_out_path.replace('.fna', '.mmseqs'), tmpdir=tmpdir
    )
    print(cmd1)
    subprocess.run(cmd1.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cluster_flanks_rep_path = flanks_out_path.replace('.fna', '.mmseqs_rep_seq.fasta')
    cluster_flanks_all_path = flanks_out_path.replace('.fna', '.mmseqs_all_seqs.fasta')
    clusters_flanks_path = flanks_out_path.replace('.fna', '.mmseqs_cluster.tsv')

    # Filter to a certain maximum number of sequences
    flank_seqs = []
    try:
        for rec in SeqIO.parse(cluster_flanks_rep_path, 'fasta'):
            flanks = recid_to_flanks[int(rec.id)]
            flank_seqs.append((int(rec.id), flanks_to_pident[flanks], flanks[0], flanks[1]))
    except FileNotFoundError:
        pass

    if random_seq:
        random.shuffle(flank_seqs)
        flank_seqs = flank_seqs[:max_seqs]
    else:
        flank_seqs = list(reversed(sorted(flank_seqs, key=lambda x: x[1])))
        flank_seqs = flank_seqs[:max_seqs]

    try:
        flanks_df = pd.DataFrame(flank_seqs)
        flanks_df.columns = ['recid', 'tpase_pident', 'left_flankseq', 'right_flankseq']
    except:
        flanks_df = pd.DataFrame(columns=['recid', 'tpase_pident', 'left_flankseq', 'right_flankseq'])

    flanks_df = flanks_df.sort_values(['tpase_pident'], ascending=False)
    flanks_df = flanks_df[['recid', 'tpase_pident', 'left_flankseq', 'right_flankseq']]
    flanks_df.to_csv(os.path.join(outdir, prefix + '.filtered.tsv'), sep='\t', index=False)

    shutil.rmtree(tmpdir)


def reverse_linfold(linfold):
    revfold = ''
    for c in linfold[::-1]:
        if c == '(':
            revfold += ')'
        elif c == ')':
            revfold += '('
        else:
            revfold += c
    return revfold


def _fold_align(filtered_tsv, aligner, max_gaps, fiveprime_user_seq, threeprime_user_seq, outdir, prefix, threads,
                revcomp_fold, keep_user_seq_cols, force):

    print("Running fold-align...")
    misc.print_variables_dictionary(locals())

    os.makedirs(outdir, exist_ok=True)
    tmpdir = os.path.join(outdir, 'tmp_' + str(random.randint(-sys.maxsize, sys.maxsize)))
    os.makedirs(tmpdir, exist_ok=True)

    filtered_df = pd.read_csv(filtered_tsv, sep='\t')
    if filtered_df is None or len(filtered_df) == 0 or filtered_df.shape[0] == 0:
        print("Filtered table is empty, exiting...")
        outdf = pd.DataFrame(columns=['recid', 'tpase_pident', 'left_flankseq', 'right_flankseq', 'left_linfold_mfe',
                                      'right_linfold_mfe', 'left_linfold', 'right_linfold', 'left_linfold_2struct',
                                      'right_linfold_2struct', 'left_flankseq_aln', 'right_flankseq_aln',
                                      'left_linfold_aln', 'right_linfold_aln', 'left_linfold_2struct_aln',
                                      'right_linfold_2struct_aln', 'left_flankseq_aln_nogaps',
                                      'right_flankseq_aln_nogaps', 'left_linfold_aln_nogaps',
                                      'right_linfold_aln_nogaps', 'left_linfold_2struct_aln_nogaps',
                                      'right_linfold_2struct_aln_nogaps'])
        outdf.to_csv(os.path.join(outdir, prefix + '.fold_align.tsv'), sep='\t', index=None)
        return

    filtered_df['recid'] = filtered_df['recid'].astype('str')


    if fiveprime_user_seq is not None:
        out = pd.DataFrame([['user', 100.000, fiveprime_user_seq, threeprime_user_seq]])
        out.columns = ['recid', 'tpase_pident', 'left_flankseq', 'right_flankseq']
        filtered_df = pd.concat([out, filtered_df])

    # Run Linearfold
    print("Running linearfold on all subflanks...")
    with Pool(threads) as pool:
        if revcomp_fold:
            left_revcomp_seqs = [reverse_complement(seq) for seq in list(filtered_df.left_flankseq)]
            right_revcomp_seqs = [reverse_complement(seq) for seq in list(filtered_df.right_flankseq)]
            left_linfold_results = pool.map(run_linearfold, left_revcomp_seqs)
            left_linfold_results = [[reverse_complement(l[0]), reverse_linfold(l[1]), l[2]] for l in
                                    left_linfold_results]
            right_linfold_results = pool.map(run_linearfold, right_revcomp_seqs)
            right_linfold_results = [[reverse_complement(l[0]), reverse_linfold(l[1]), l[2]] for l in
                                     right_linfold_results]
        else:
            left_linfold_results = pool.map(run_linearfold, list(filtered_df.left_flankseq))
            right_linfold_results = pool.map(run_linearfold, list(filtered_df.right_flankseq))

    filtered_df['left_linfold_mfe'] = [l[2] for l in left_linfold_results]
    filtered_df['right_linfold_mfe'] = [l[2] for l in right_linfold_results]
    filtered_df['left_linfold'] = [l[1] for l in left_linfold_results]
    filtered_df['right_linfold'] = [l[1] for l in right_linfold_results]

    left_bulge_struc = []
    right_bulge_struc = []
    for left_linfold, right_linfold in zip(filtered_df.left_linfold, filtered_df.right_linfold):
        left_bg = bulge_graph.BulgeGraph.from_dotbracket(left_linfold)
        right_bg = bulge_graph.BulgeGraph.from_dotbracket(right_linfold)

        left_secondary_struct = "".join([left_bg.get_elem(idx + 1)[0] for idx in range(len(left_linfold))]).upper()
        right_secondary_struct = "".join([right_bg.get_elem(idx + 1)[0] for idx in range(len(right_linfold))]).upper()

        left_bulge_struc.append(left_secondary_struct)
        right_bulge_struc.append(right_secondary_struct)

    filtered_df['left_linfold_2struct'] = left_bulge_struc
    filtered_df['right_linfold_2struct'] = right_bulge_struc

    lf_filtered_path = os.path.join(tmpdir, 'lf.filt.fna')
    rf_filtered_path = os.path.join(tmpdir, 'rf.filt.fna')

    num_lf_seqs = 0
    outfile = open(lf_filtered_path, 'w')
    for recid, flankseq in zip(filtered_df.recid, filtered_df.left_flankseq):
        num_lf_seqs += 1
        print(">" + str(recid), flankseq, sep='\n', file=outfile)
    outfile.close()

    num_rf_seqs = 0
    outfile = open(rf_filtered_path, 'w')
    for recid, flankseq in zip(filtered_df.recid, filtered_df.right_flankseq):
        num_rf_seqs += 1
        print(">" + str(recid), flankseq, sep='\n', file=outfile)
    outfile.close()

    # Run MUSCLE
    muscle_lf_align_out = os.path.join(tmpdir, 'lf.filt.aln.fna')
    muscle_rf_align_out = os.path.join(tmpdir, 'rf.filt.aln.fna')

    # Run MUSCLE
    if aligner == 'muscle':
        print("Aligning with MUSCLE...")

        if (num_lf_seqs < 100):
            cmd = '{muscle} -align {fna} -output {out}'.format(muscle=MUSCLE_PATH, fna=lf_filtered_path,
                                                               out=muscle_lf_align_out)
        else:
            cmd = '{muscle} -super5 {fna} -output {out}'.format(muscle=MUSCLE_PATH, fna=lf_filtered_path,
                                                                out=muscle_lf_align_out)
        print(cmd)
        subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        if (num_rf_seqs < 100):
            cmd = '{muscle} -align {fna} -output {out}'.format(muscle=MUSCLE_PATH, fna=rf_filtered_path,
                                                               out=muscle_rf_align_out)
        else:
            cmd = '{muscle} -super5 {fna} -output {out}'.format(muscle=MUSCLE_PATH, fna=rf_filtered_path,
                                                                out=muscle_rf_align_out)
        print(cmd)
        subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    elif aligner == 'mafft-qinsi':
        print("Aligning with MAFFT-QINSI...")
        muscle_reticle_align_out = os.path.join(tmpdir, 'reticles.filt.aln.fna')
        cmd = '{mafft} --thread {t} --maxiterate 1000 {fna} > {out}'.format(mafft=MAFFT_QINSI_PATH, t=threads,
                                                                            fna=lf_filtered_path,
                                                                            out=muscle_lf_align_out)
        subprocess.run(cmd, shell=True, capture_output=True, check=True)
        cmd = '{mafft} --thread {t} --maxiterate 1000 {fna} > {out}'.format(
            mafft=MAFFT_QINSI_PATH, t=threads, fna=rf_filtered_path, out=muscle_rf_align_out
        )
        subprocess.run(cmd, shell=True, capture_output=True, check=True)

    elif aligner == 'mafft-xinsi':
        print("Aligning with MAFFT-XINSI...")
        muscle_reticle_align_out = os.path.join(tmpdir, 'reticles.filt.aln.fna')
        cmd = '{mafft} --thread {t} --maxiterate 1000 {fna} > {out}'.format(mafft=MAFFT_XINSI_PATH, t=threads,
                                                                            fna=lf_filtered_path,
                                                                            out=muscle_lf_align_out)
        subprocess.run(cmd, shell=True, capture_output=True, check=True)  # , capture_output=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cmd = '{mafft} --thread {t} --maxiterate 1000 {fna} > {out}'.format(mafft=MAFFT_XINSI_PATH, t=threads,
                                                                            fna=rf_filtered_path,
                                                                            out=muscle_rf_align_out)
        subprocess.run(cmd, shell=True, capture_output=True, check=True)  # , capture_output=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Make data frame
    lf_recid_to_aln = dict()
    try:
        for rec in SeqIO.parse(muscle_lf_align_out, 'fasta'):
            lf_recid_to_aln[rec.id] = str(rec.seq)
    except FileNotFoundError:
        pass

    rf_recid_to_aln = dict()
    try:
        for rec in SeqIO.parse(muscle_rf_align_out, 'fasta'):
            rf_recid_to_aln[rec.id] = str(rec.seq)
    except FileNotFoundError:
        pass

    filtered_df['left_flankseq_aln'] = [lf_recid_to_aln[recid] for recid in filtered_df.recid]
    filtered_df['right_flankseq_aln'] = [rf_recid_to_aln[recid] for recid in filtered_df.recid]

    filtered_df['left_linfold_aln'] = [map_linfold_to_aln(aln, linfold) for aln, linfold in
                                       zip(filtered_df.left_flankseq_aln, filtered_df.left_linfold)]
    filtered_df['right_linfold_aln'] = [map_linfold_to_aln(aln, linfold) for aln, linfold in
                                        zip(filtered_df.right_flankseq_aln, filtered_df.right_linfold)]

    filtered_df['left_linfold_2struct_aln'] = [map_linfold_to_aln(aln, linfold) for aln, linfold in
                                               zip(filtered_df.left_flankseq_aln, filtered_df.left_linfold_2struct)]
    filtered_df['right_linfold_2struct_aln'] = [map_linfold_to_aln(aln, linfold) for aln, linfold in
                                                zip(filtered_df.right_flankseq_aln, filtered_df.right_linfold_2struct)]

    left_gap_counts = defaultdict(int)
    total_left = 0
    for aln in filtered_df.left_flankseq_aln:
        total_left += 1
        for i, c in enumerate(aln):
            if c == '-':
                left_gap_counts[i] += 1

    left_gappy_positions = set()
    for i in left_gap_counts:
        if left_gap_counts[i] / total_left * 100 > max_gaps:
            left_gappy_positions.add(i)

    if keep_user_seq_cols:
        user_seq_aln = list(filtered_df[filtered_df.recid == 'user'].left_flankseq_aln)[0]
        for i, c in enumerate(user_seq_aln):
            if c == '-':
                continue
            left_gappy_positions.discard(i)

    right_gap_counts = defaultdict(int)
    total_right = 0
    for aln in filtered_df.right_flankseq_aln:
        total_right += 1
        for i, c in enumerate(aln):
            if c == '-':
                right_gap_counts[i] += 1

    right_gappy_positions = set()
    for i in right_gap_counts:
        if right_gap_counts[i] / total_right * 100 > max_gaps:
            right_gappy_positions.add(i)

    if keep_user_seq_cols:
        user_seq_aln = list(filtered_df[filtered_df.recid == 'user'].right_flankseq_aln)[0]
        for i, c in enumerate(user_seq_aln):
            if c == '-':
                continue
            right_gappy_positions.discard(i)

    filtered_df['left_flankseq_aln_nogaps'] = [remove_gaps(aln, left_gappy_positions) for aln in
                                               filtered_df.left_flankseq_aln]
    filtered_df['right_flankseq_aln_nogaps'] = [remove_gaps(aln, right_gappy_positions) for aln in
                                                filtered_df.right_flankseq_aln]
    filtered_df['left_linfold_aln_nogaps'] = [remove_gaps(aln, left_gappy_positions) for aln in
                                              filtered_df.left_linfold_aln]
    filtered_df['right_linfold_aln_nogaps'] = [remove_gaps(aln, right_gappy_positions) for aln in
                                               filtered_df.right_linfold_aln]
    filtered_df['left_linfold_2struct_aln_nogaps'] = [remove_gaps(aln, left_gappy_positions) for aln in
                                                      filtered_df.left_linfold_2struct_aln]
    filtered_df['right_linfold_2struct_aln_nogaps'] = [remove_gaps(aln, right_gappy_positions) for aln in
                                                       filtered_df.right_linfold_2struct_aln]
    filtered_df.to_csv(os.path.join(outdir, prefix + '.fold_align.tsv'), sep='\t', index=None)

    shutil.rmtree(tmpdir)


def get_structure_perc(structures, dot):
    struc_counts = defaultdict(lambda: defaultdict(int))
    for seq, dot in zip(structures, dot):
        for i, c in enumerate(seq):
            if c == 'S':
                if dot[i] == '(':
                    struc_counts[i]['5S'] += 1
                elif dot[i] == ')':
                    struc_counts[i]['3S'] += 1
                else:
                    print("STEM ERROR")
                    sys.exit()
            else:
                struc_counts[i][c] += 1

    for i in struc_counts:
        for c in struc_counts[i]:
            struc_counts[i][c] = struc_counts[i][c] / len(structures) * 100

    return struc_counts


def get_structure_regions(struc_perc, struc, cutoff):
    consec_structures = []
    this_consec = []
    for i in range(len(struc_perc)):
        perc = struc_perc[i][struc]
        if perc > 50:
            if len(this_consec) == 0 or this_consec[-1] == i - 1:
                this_consec.append(i)
            else:
                consec_structures.append(list(this_consec))
                this_consec = [i]
    if len(this_consec) > 0:
        consec_structures.append(this_consec)

    return consec_structures


def calc_stem_statistics(hairpin_regions, five_stem_regions, three_stem_regions):
    stats = []
    for hp in hairpin_regions:
        left_five_stem_regions = []
        right_five_stem_regions = []
        for stm in five_stem_regions:
            if stm[0] < hp[0] and stm[-1] < hp[-1]:
                left_five_stem_regions.append(stm)
            elif stm[-1] > hp[-1] and stm[0] > hp[0]:
                right_five_stem_regions.append(stm)

        left_three_stem_regions = []
        right_three_stem_regions = []
        for stm in three_stem_regions:
            if stm[0] < hp[0] and stm[-1] < hp[-1]:
                left_three_stem_regions.append(stm)
            elif stm[-1] > hp[-1] and stm[0] > hp[0]:
                right_three_stem_regions.append(stm)

        left_stems_clean = []
        if len(left_three_stem_regions) == 0:
            left_stems_clean = list(left_five_stem_regions)
        else:
            for stm in left_five_stem_regions:
                if stm[0] > left_three_stem_regions[-1][-1]:
                    left_stems_clean.append(stm)

        right_stems_clean = []
        if len(right_five_stem_regions) == 0:
            right_stems_clean = list(right_three_stem_regions)
        else:
            for stm in right_three_stem_regions:
                if stm[-1] < right_five_stem_regions[0][0]:
                    right_stems_clean.append(stm)

        left_stem_total = sum([1 for e in left_stems_clean for j in e])
        right_stem_total = sum([1 for e in right_stems_clean for j in e])

        try:
            left_seq_length = hp[0] - left_stems_clean[0][0]
        except IndexError:
            left_seq_length = 0

        try:
            right_seq_length = right_stems_clean[-1][-1] - hp[-1]
        except:
            right_seq_length = 0

        left_interior_total = 0
        left_interior_max_consec = 0
        if len(left_stems_clean) > 0:
            for i in range(1, len(left_stems_clean)):
                inter_len = left_stems_clean[i][0] - left_stems_clean[i - 1][-1]
                left_interior_total += inter_len
                if inter_len > left_interior_max_consec:
                    left_interior_max_consec = inter_len

        right_interior_total = 0
        right_interior_max_consec = 0
        if len(right_stems_clean) > 0:
            for i in range(1, len(right_stems_clean)):
                inter_len = right_stems_clean[i][0] - right_stems_clean[i - 1][-1]
                right_interior_total += inter_len
                if inter_len > right_interior_max_consec:
                    right_interior_max_consec = inter_len

        try:
            five_hairpin_stem_length = len(left_stems_clean[-1])
        except IndexError:
            five_hairpin_stem_length = 0

        try:
            three_hairpin_stem_length = len(right_stems_clean[0])
        except IndexError:
            three_hairpin_stem_length = 0

        try:
            five_stem_start = left_stems_clean[0][0]
        except:
            five_stem_start = hp[0]
        try:
            three_stem_end = right_stems_clean[-1][-1] + 1
        except:
            three_stem_end = hp[-1]

        stats.append([hp[0], hp[-1], max(hp) - min(hp) + 1, five_hairpin_stem_length, three_hairpin_stem_length,
                      five_stem_start, three_stem_end, left_stem_total, right_stem_total,
                      left_seq_length, right_seq_length, left_interior_total, right_interior_total,
                      left_interior_max_consec, right_interior_max_consec])

    stats = pd.DataFrame(stats)
    stats.columns = ['hairpin_start', 'hairpin_end', 'hairpin_length', 'five_hairpin_stem_length',
                     'three_hairpin_stem_length', 'five_stem_start', 'three_stem_end', 'left_stem_total',
                     'right_stem_total', 'left_seq_length', 'right_seq_length', 'left_interior_total',
                     'right_interior_total',
                     'left_interior_max_consec', 'right_interior_max_consec']

    return stats


def get_boundaries(midpoint, struc_perc, stem_start, stem_end, structure_cutoff):
    print('midpoint:', midpoint)
    print('stem_start:', stem_start)
    print('stem_end:', stem_end)
    print('structure_cutoff:', structure_cutoff)

    five_stem_regions = get_structure_regions(struc_perc, '5S', structure_cutoff)
    three_stem_regions = get_structure_regions(struc_perc, '3S', structure_cutoff)

    five_stem_regions = list(
        reversed([j for e in five_stem_regions for j in e if e[-1] < midpoint and e[0] >= stem_start]))
    three_stem_regions = [j for e in three_stem_regions for j in e if e[0] > midpoint and e[-1] < stem_end]

    num_stems = min([len(five_stem_regions), len(three_stem_regions)])
    if num_stems > 0:
        indices = set()
        for i in range(num_stems):
            indices.add(five_stem_regions[i])
            indices.add(three_stem_regions[i])
        indices = list(indices)
        bound_start, bound_end = min(indices), max(indices) + 1

    else:
        bound_start, bound_end = stem_start, stem_end

    return bound_start, bound_end, num_stems


def extract_consensus_structure(midpoint, start, end, flank, struc_perc, structure_cutoff):
    five_stem_regions = get_structure_regions(struc_perc, '5S', structure_cutoff)
    three_stem_regions = get_structure_regions(struc_perc, '3S', structure_cutoff)

    five_stem_regions = [j for e in five_stem_regions for j in e if e[-1] < midpoint]
    three_stem_regions = [j for e in three_stem_regions for j in e if e[0] > midpoint]

    five_stem_regions = [i for i in five_stem_regions if i >= start]
    three_stem_regions = [i for i in three_stem_regions if i < end]

    flank_start = start - flank
    if flank_start < 0:
        flank_start = 0

    flank_end = end + flank
    if flank_end > max(struc_perc.keys()):
        flank_end = max(struc_perc.keys())

    consensus = ""
    for i in range(flank_start, flank_end):
        if i in five_stem_regions:
            consensus += "<"
        elif i in three_stem_regions:
            consensus += ">"
        else:
            consensus += "."

    return consensus, flank_start, flank_end


def extract_alignments(aln, start, end):
    out = []
    for s in aln:
        out.append(s[start:end])

    return out


def write_stockholm_alignment(aln, consensus, outpath):
    cons_chunks = []
    for j in range(0, len(consensus), 50):
        cons_chunks.append(consensus[j:j + 50])

    aln_chunks = []
    for s in aln:
        chunks = []
        for j in range(0, len(s), 50):
            chunks.append(s[j:j + 50])
        aln_chunks.append(chunks)

    outfile = open(outpath, 'w')
    print("# STOCKHOLM 1.0\n", file=outfile)
    for i in range(len(cons_chunks)):
        for j in range(len(aln_chunks)):
            name = 'seq' + str(j + 1)
            print(name + " " * (24 - len(name)) + aln_chunks[j][i], file=outfile)

        if i != len(cons_chunks) - 1:
            name = "#=GC SS_cons"
            print(name + " " * (24 - len(name)) + cons_chunks[i] + "\n", file=outfile)
        else:
            name = "#=GC SS_cons"
            print(name + " " * (24 - len(name)) + cons_chunks[i] + "\n//", file=outfile)

    outfile.close()


def _detect(foldalign, max_seq_gaps, structure_cutoff, drna_flank_length, outdir, prefix, force):
    os.makedirs(outdir, exist_ok=True)

    foldalign_df = pd.read_csv(foldalign, sep='\t')

    keep_indices = set()
    remove_indices = set()
    for index, seq in enumerate(foldalign_df.linfold_2struct_aln_nogaps):

        missing_count = seq.count("-")
        if missing_count / len(seq) * 100 > max_seq_gaps:
            remove_indices.add(index)
        else:
            keep_indices.add(index)

    print("Removing %d sequences due to high gappyness..." % len(remove_indices))
    foldalign_df = foldalign_df.iloc[list(keep_indices)]
    foldalign_left = foldalign_df[foldalign_df.flankside == 'left']
    foldalign_right = foldalign_df[foldalign_df.flankside == 'right']

    # do left
    left_struc_perc = get_structure_perc(foldalign_left.linfold_2struct_aln_nogaps, foldalign_left.linfold_aln_nogaps)
    hairpin_regions = get_structure_regions(left_struc_perc, 'H', structure_cutoff)

    print("Detected %d candidate hairpin regions in the left CDS flank..." % len(hairpin_regions))
    if len(hairpin_regions) > 0:
        five_stem_regions = get_structure_regions(left_struc_perc, '5S', structure_cutoff)
        three_stem_regions = get_structure_regions(left_struc_perc, '3S', structure_cutoff)

        left_stem_stats = calc_stem_statistics(hairpin_regions, five_stem_regions, three_stem_regions)
    else:
        left_stem_stats = pd.DataFrame(columns=[
            'hairpin_start', 'hairpin_end', 'hairpin_length', 'five_hairpin_stem_length',
            'three_hairpin_stem_length', 'five_stem_start', 'three_stem_end', 'left_stem_total',
            'right_stem_total', 'left_seq_length', 'right_seq_length', 'left_interior_total', 'right_interior_total',
            'left_interior_max_consec', 'right_interior_max_consec'])

    # do right
    right_struc_perc = get_structure_perc(foldalign_right.linfold_2struct_aln_nogaps,
                                          foldalign_right.linfold_aln_nogaps)
    hairpin_regions = get_structure_regions(right_struc_perc, 'H', structure_cutoff)

    print("Detected %d candidate hairpin regions in the right CDS flank..." % len(hairpin_regions))
    if len(hairpin_regions) > 0:
        five_stem_regions = get_structure_regions(right_struc_perc, '5S', structure_cutoff)
        three_stem_regions = get_structure_regions(right_struc_perc, '3S', structure_cutoff)

        right_stem_stats = calc_stem_statistics(hairpin_regions, five_stem_regions, three_stem_regions)
    else:
        right_stem_stats = pd.DataFrame(columns=[
            'hairpin_start', 'hairpin_end', 'hairpin_length', 'five_hairpin_stem_length',
            'three_hairpin_stem_length', 'five_stem_start', 'three_stem_end', 'left_stem_total',
            'right_stem_total', 'left_seq_length', 'right_seq_length', 'left_interior_total', 'right_interior_total',
            'left_interior_max_consec', 'right_interior_max_consec'])

    left_stem_stats['flankside'] = 'left'
    right_stem_stats['flankside'] = 'right'

    stem_stats = pd.concat([left_stem_stats, right_stem_stats])
    stem_stats['score'] = np.minimum(np.abs(5 - stem_stats.hairpin_length), 10) + np.minimum(
        np.abs(30 - stem_stats.left_seq_length), 10) + np.minimum(np.abs(30 - stem_stats.right_seq_length),
                                                                  10) + np.minimum(
        np.abs(12 - stem_stats.left_interior_max_consec), 10) + np.minimum(
        np.abs(12 - stem_stats.right_interior_max_consec), 10)

    stem_stats = stem_stats.sort_values('score').reset_index()

    if stem_stats.shape[0] == 0:
        out = open(os.path.join(outdir, prefix + '.detect.info.tsv'), 'w')
        out.close()

        out = open(os.path.join(outdir, prefix + '.detect.aln.sto'), 'w')
        out.close()

        return False

    # Pick top candidate
    top_position = stem_stats.iloc[0]

    if top_position.flankside == 'left':
        print("Top candidate found in left flank...")
        struc_perc = left_struc_perc
        nuc_aln = foldalign_left.flankseq_aln_nogaps
    else:
        print("Top candidate found in right flank...")
        struc_perc = right_struc_perc
        nuc_aln = foldalign_right.flankseq_aln_nogaps

    # Find boundaries
    drna_start, drna_end, num_stems = get_boundaries(int((top_position.hairpin_start + top_position.hairpin_end) / 2),
                                                     struc_perc, top_position.five_stem_start,
                                                     top_position.three_stem_end, structure_cutoff)

    # Get dRNA structure

    drna_consensus, drna_extract_start, drna_extract_end = extract_consensus_structure(
        int((top_position.hairpin_start + top_position.hairpin_end) / 2), drna_start, drna_end, drna_flank_length,
        struc_perc, structure_cutoff)

    extracted_alignments = extract_alignments(nuc_aln, drna_extract_start, drna_extract_end)
    write_stockholm_alignment(extracted_alignments, drna_consensus, os.path.join(outdir, prefix + '.detect.aln.sto'))

    outfile = open(os.path.join(outdir, prefix + '.detect.info.tsv'), 'w')
    print("field", "value", sep='\t', file=outfile)
    print("flankside", top_position.flankside, sep='\t', file=outfile)
    print("score", top_position.score, sep='\t', file=outfile)
    print("hairpin_start", top_position.hairpin_start, sep='\t', file=outfile)
    print("hairpin_end", top_position.hairpin_end + 1, sep='\t', file=outfile)
    print("struc_stem_start", drna_start, sep='\t', file=outfile)
    print("struc_stem_end", drna_end, sep='\t', file=outfile)
    print("struc_flank_start", drna_extract_start, sep='\t', file=outfile)
    print("struc_flank_end", drna_extract_end, sep='\t', file=outfile)
    print("num_stems", num_stems, sep='\t', file=outfile)
    print("five_hairpin_stem_length", top_position.five_hairpin_stem_length, sep='\t', file=outfile)
    print("three_hairpin_stem_length", top_position.three_hairpin_stem_length, sep='\t', file=outfile)
    print("left_stem_total", top_position.left_stem_total, sep='\t', file=outfile)
    print("right_stem_total", top_position.right_stem_total, sep='\t', file=outfile)
    print("left_seq_length", top_position.left_seq_length, sep='\t', file=outfile)
    print("right_seq_length", top_position.right_seq_length, sep='\t', file=outfile)
    print("left_interior_total", top_position.left_interior_total, sep='\t', file=outfile)
    print("right_interior_total", top_position.right_interior_total, sep='\t', file=outfile)
    print("left_interior_max_consec", top_position.left_interior_max_consec, sep='\t', file=outfile)
    print("right_interior_max_consec", top_position.right_interior_max_consec, sep='\t', file=outfile)
    outfile.close()


def extract_aln_range(aln, unalign_seq_len):
    start = -1
    end = -1
    unalign_count = 0
    for i, c in enumerate(aln):
        if c != '-' and start == -1:
            start = i
            unalign_count += 1
        elif c != '-' and unalign_count == unalign_seq_len - 1:
            end = i + 1
            break
        elif c != '-':
            unalign_count += 1
    return start, end


def _reticle_refold(foldalign_tsv, aligner, max_gaps, tic_length, outdir, prefix, threads, keep_user_seq_cols,
                    revcomp_fold, force):
    os.makedirs(outdir, exist_ok=True)
    tmpdir = os.path.join(outdir, 'tmp_' + str(random.randint(-sys.maxsize, sys.maxsize)))
    os.makedirs(tmpdir, exist_ok=True)

    foldalign_df = pd.read_csv(foldalign_tsv, sep='\t')

    left_user_seq = list(foldalign_df[foldalign_df.recid == 'user'].left_flankseq)[0]
    right_user_seq = list(foldalign_df[foldalign_df.recid == 'user'].right_flankseq)[0]

    left_user_seq_len = len(left_user_seq)
    right_user_seq_len = len(right_user_seq) - tic_length

    left_user_seq_aln = list(foldalign_df[foldalign_df.recid == 'user'].left_flankseq_aln)[0]
    right_user_seq_aln = list(foldalign_df[foldalign_df.recid == 'user'].right_flankseq_aln)[0]

    left_aln_range = extract_aln_range(left_user_seq_aln, left_user_seq_len)
    right_aln_range = extract_aln_range(right_user_seq_aln, right_user_seq_len)

    reticles = []
    for recid, tpase_pident, right_aln, left_aln in zip(foldalign_df.recid, foldalign_df.tpase_pident,
                                                        foldalign_df.right_flankseq_aln,
                                                        foldalign_df.left_flankseq_aln):
        reticles.append((recid, tpase_pident, (right_aln[right_aln_range[0]:right_aln_range[1]] + left_aln[
                                                                                                  left_aln_range[0]:
                                                                                                  left_aln_range[
                                                                                                      1]]).replace("-",
                                                                                                                   "")))

    reticles_df = pd.DataFrame(reticles)
    reticles_df.columns = ['recid', 'tpase_pident', 'reticle_seq']

    # Run Linearfold
    print("Running linearfold on all reticles...")
    with Pool(threads) as pool:
        if revcomp_fold:
            reticle_revcomp_seqs = [reverse_complement(seq) for seq in list(reticles_df.reticle_seq)]
            reticle_linfold_results = pool.map(run_linearfold, reticle_revcomp_seqs)
            reticle_linfold_results = [[reverse_complement(l[0]), reverse_linfold(l[1]), l[2]] for l in
                                       reticle_linfold_results]
        else:
            reticle_linfold_results = pool.map(run_linearfold, list(reticles_df.reticle_seq))

    reticles_df['reticle_linfold_mfe'] = [l[2] for l in reticle_linfold_results]
    reticles_df['reticle_linfold'] = [l[1] for l in reticle_linfold_results]

    reticle_bulge_struc = []
    for linfold in reticles_df.reticle_linfold:
        reticle_bg = bulge_graph.BulgeGraph.from_dotbracket(linfold)

        reticle_secondary_struct = "".join([reticle_bg.get_elem(idx + 1)[0] for idx in range(len(linfold))]).upper()

        reticle_bulge_struc.append(reticle_secondary_struct)

    reticles_df['reticle_linfold_2struct'] = reticle_bulge_struc

    reticle_filtered_path = os.path.join(tmpdir, 'reticles.filt.fna')

    num_reticle_seqs = 0
    outfile = open(reticle_filtered_path, 'w')
    for recid, retseq in zip(reticles_df.recid, reticles_df.reticle_seq):
        print(">" + str(recid), retseq, sep='\n', file=outfile)
    outfile.close()

    # Run MUSCLE
    if aligner == 'muscle':
        print("Aligning with MUSCLE...")
        muscle_reticle_align_out = os.path.join(tmpdir, 'reticles.filt.aln.fna')

        if (num_reticle_seqs < 100):
            cmd = '{muscle} -align {fna} -output {out}'.format(muscle=MUSCLE_PATH, fna=reticle_filtered_path,
                                                               out=muscle_reticle_align_out)
        else:
            cmd = '{muscle} -super5 {fna} -output {out}'.format(muscle=MUSCLE_PATH, fna=reticle_filtered_path,
                                                                out=muscle_reticle_align_out)
        print(cmd)
        subprocess.run(cmd.split(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    elif aligner == 'mafft-qinsi':
        print("Aligning with MAFFT-QINSI...")
        muscle_reticle_align_out = os.path.join(tmpdir, 'reticles.filt.aln.fna')
        cmd = '{mafft} --maxiterate 1000 {fna} > {out}'.format(mafft=MAFFT_QINSI_PATH, fna=reticle_filtered_path,
                                                               out=muscle_reticle_align_out)
        subprocess.run(cmd, shell=True, capture_output=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Make data frame
    reticle_recid_to_aln = dict()
    try:
        for rec in SeqIO.parse(muscle_reticle_align_out, 'fasta'):
            reticle_recid_to_aln[rec.id] = str(rec.seq)
    except FileNotFoundError:
        pass

    reticles_df['reticle_aln'] = [reticle_recid_to_aln[recid] for recid in reticles_df.recid]

    reticles_df['reticle_linfold_aln'] = [map_linfold_to_aln(aln, linfold) for aln, linfold in
                                          zip(reticles_df.reticle_aln, reticles_df.reticle_linfold)]
    reticles_df['reticle_linfold_2struct_aln'] = [map_linfold_to_aln(aln, linfold) for aln, linfold in
                                                  zip(reticles_df.reticle_aln, reticles_df.reticle_linfold_2struct)]

    reticle_gap_counts = defaultdict(int)
    total_reticle = 0
    for aln in reticles_df.reticle_aln:
        total_reticle += 1
        for i, c in enumerate(aln):
            if c == '-':
                reticle_gap_counts[i] += 1

    reticle_gappy_positions = set()
    for i in reticle_gap_counts:
        if reticle_gap_counts[i] / total_reticle * 100 > max_gaps:
            reticle_gappy_positions.add(i)

    if keep_user_seq_cols:
        user_reticle_aln = list(reticles_df[reticles_df.recid == 'user'].reticle_aln)[0]
        for i, c in enumerate(user_reticle_aln):
            if c == '-':
                continue
            reticle_gappy_positions.discard(i)

    reticles_df['reticle_aln_nogaps'] = [remove_gaps(aln, reticle_gappy_positions) for aln in reticles_df.reticle_aln]
    reticles_df['reticle_linfold_aln_nogaps'] = [remove_gaps(aln, reticle_gappy_positions) for aln in
                                                 reticles_df.reticle_linfold_aln]
    reticles_df['reticle_linfold_2struct_aln_nogaps'] = [remove_gaps(aln, reticle_gappy_positions) for aln in
                                                         reticles_df.reticle_linfold_2struct_aln]
    reticles_df.to_csv(os.path.join(outdir, prefix + '.reticle_refold.tsv'), sep='\t', index=None)

    shutil.rmtree(tmpdir)


def _retrieve_elements(blastp_tsv, min_pident, min_alnlen, max_unique_proteins, outdir, prefix, force):
    os.makedirs(outdir, exist_ok=True)

    keep_loci = set()

    blastp = pd.read_csv(blastp_tsv, sep='\t')
    blastp = blastp[blastp.pident > min_pident]
    blastp = blastp[(blastp.alnlen / blastp.qlen) * 100 > min_alnlen]
    blastp = blastp[(blastp.alnlen / blastp.plen) * 100 > min_alnlen]
    blastp = blastp.sort_values('pident', ascending=False)

    keep_proteins = set()
    for prot in blastp.p100_id:
        if len(keep_proteins) >= max_unique_proteins:
            break
        keep_proteins.add(prot)

    if len(keep_proteins) == 0:
        print("No proteins met filtering thresholds.")
        return False
    else:
        print("{} proteins met filtering thresholds...".format(len(keep_proteins)))

    p100_to_pident = dict()
    already = set()
    for p100, pident in zip(blastp.p100_id, blastp.pident):
        if p100 in already:
            continue
        already.add(p100)
        p100_to_pident[p100] = pident

    # Get elements

    outfile = open(os.path.join(outdir, prefix + '.retrieve_elements.tsv'), 'w')
    print("locus_id", "element_id", "p100_id", "search_iteration", "pident", "core", "donor_seq", "target_seq",
          "elemseq", sep='\t', file=outfile)
    with sqlite3.connect(IS110_DB) as conn:

        cur = conn.cursor()

        for prot in tqdm(keep_proteins, desc='Protein flanks retrieved'):

            cmd = f'SELECT locus_id, locus_seq FROM loci WHERE p100_id="{prot}"'
            prot_seq = None
            locus_ids = set()
            for res in cur.execute(cmd):
                locus_ids.add(res)

            if len(locus_ids) == 0:
                continue

            for lid, lseq in locus_ids:
                cmd = f'SELECT * FROM loci_annot WHERE locus_id="{lid}";'

                elements = defaultdict(lambda: {'cl': None, 'ec': None})
                for res in cur.execute(cmd):

                    elemhash = res[-1].split('=')[-1]
                    if res[1] == 'LTIM':
                        corelength = res[3] - res[2]
                        elements[elemhash]['cl'] = corelength
                    if res[1] == 'IS110':
                        elements[elemhash]['ec'] = res[2], res[3]

                for e in elements:
                    cmd = f'SELECT search_iteration FROM elements WHERE element_id="{e}"'
                    for res in cur.execute(cmd):
                        siter = res[0]
                        break
                    info = elements[e]
                    core = lseq[info['ec'][0] - info['cl']:info['ec'][0]]
                    elemseq = lseq[info['ec'][0]:info['ec'][1]]
                    donor_seq = lseq[info['ec'][1] - 20:info['ec'][1]] + core + lseq[info['ec'][0]:info['ec'][0] + 20]
                    target_seq = lseq[info['ec'][0] - (20 + len(core)):info['ec'][0] - len(core)] + core + lseq[
                                                                                                           info['ec'][
                                                                                                               1] + len(
                                                                                                               core):
                                                                                                           info['ec'][
                                                                                                               1] + 20 + len(
                                                                                                               core)]
                    print(lid, e, prot, siter, p100_to_pident[prot], core, donor_seq, target_seq, elemseq, sep='\t',
                          file=outfile)
    outfile.close()
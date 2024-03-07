import subprocess, os, re, shutil
from bridgerna2024 import *
from bridgerna2024 import sctools, misc
import pysam
from collections import defaultdict
from Bio.Seq import reverse_complement
from Bio import SeqIO
import pandas as pd
import tempfile
from hashlib import md5 as hashlib_md5
from base64 import b64encode
import random
import gzip

def hash_string(string):
    return b64encode(hashlib_md5(string.encode('utf-8')).digest()).decode('utf-8').rstrip('=')

def _run_nanoporepipeline_snakemake(workdir, threads):
    """Run the RNA-seq snakemake pipeline."""
    subprocess.run(['snakemake', '--snakefile', NANOPOREPIPELINE_SNAKEFILE_PATH, '--configfile',
                    NANOPOREPIPELINE_SNAKEMAKE_CONFIG_PATH, '-j', str(threads),
                    '--config',  'workdir={}'.format(workdir)], check=True)

def read_length_distribution(infastq, outtsv):
    lengthcounts = defaultdict(int)
    for read in misc.fastq_reader(infastq):
        lengthcounts[len(read[1])] += 1
    outtsv = open(outtsv, 'w')
    print('read_length', 'count', sep='\t', file=outtsv)
    for length, count in reversed(sorted(lengthcounts.items(), key=lambda x: x[1])):
        print(length, count, sep='\t', file=outtsv)
    outtsv.close()

def parse_bridgerna_guides(inbam, bridgerna_fasta, outtsv):

    # Get bridgeRNA sequence
    bridgerna = None
    for rec in SeqIO.parse(bridgerna_fasta, 'fasta'):
        bridgerna = str(rec.seq)

    regions = ['ltg', 'rtg']
    position2region = dict()
    region2position = {'ltg':None, 'rtg':None}
    # Specify the different variable regions of the bridgeRNA
    for match in re.finditer(r'N+', bridgerna):
        if region2position['ltg'] is None:
            region2position['ltg'] = list(range(match.start(), match.end()))
        else:
            region2position['rtg'] = list(range(match.start(), match.end()))

    for region in region2position:
        for pos in region2position[region]:
            position2region[pos] = region

    tsv_out = open(outtsv, 'w')
    print("read_hash_id", 'ltg', 'rtg', sep='\t', file=tsv_out)

    bam = pysam.AlignmentFile(inbam, 'rb')
    for read in bam:

        # Skip unmapped reads
        if read.is_unmapped:
            continue
        # Skip secondary alignments
        if read.is_secondary:
            continue
        # Skip supplementary alignments
        if read.is_supplementary:
            continue

        read_query_seq = read.query_sequence
        read_query_quals = read.query_qualities

        read_info = {'pos2nuc': dict(), 'pos2qual': dict()}

        read_pairs = read.get_aligned_pairs()
        for query_pos, ref_pos in read_pairs:
            if query_pos is None:
                continue
            if ref_pos not in position2region:
                continue
            try:
                read_info['pos2nuc'][ref_pos] = read_query_seq[query_pos]
                read_info['pos2qual'][ref_pos] = read_query_quals[query_pos]
            except Exception as e:
                print(query_pos)
                raise e

        # Get consensus sequence for each region
        region_consenses = dict()
        for region in regions:

            region_consensus = ''
            for rpos in region2position[region]:
                if rpos not in read_info['pos2nuc']:
                    region_consensus += 'N'
                else:
                    nuc, qual = "N", 0

                    if rpos in read_info['pos2nuc']:
                        nuc, qual = read_info['pos2nuc'][rpos], read_info['pos2qual'][rpos]

                    region_consensus += nuc

            region_consenses[region] = region_consensus

        print(read.query_name, region_consenses['ltg'], region_consenses['rtg'], sep='\t', file=tsv_out)

    bam.close()
    tsv_out.close()

def count_bridgerna_guides(guides, outtsv):

    # Extract sample name from file path
    outfile = open(outtsv, 'w')
    print("sample", "ltg", "rtg", "count", sep='\t', file=outfile)
    for gf in guides:
        sample = os.path.basename(gf).replace('.guides.tsv', '')
        guide_counts = defaultdict(int)
        with open(gf) as infile:
            infile.readline()
            for line in infile:
                line = line.strip().split('\t')
                guide = tuple(line[1:])
                if "N" in guide[0] or "N" in guide[1]:
                    continue
                guide_counts[guide] += 1
        guide_counts = list(sorted(guide_counts.items(), key=lambda x: x[1], reverse=True))
        for g, c in guide_counts:
            print(sample, *g, c, sep='\t', file=outfile)
    outfile.close()

def combine_fastq(infastq, outfastq):

    outfile = gzip.open(outfastq, 'wt')
    for f in infastq:
        for name, seq, spacer, qual in misc.fastq_reader(f):
            print(name, seq, spacer, qual, sep='\n', file=outfile)
    outfile.close()

def combine_downsample_fastq(sample, infastqs, sample_info, outfastq):

    random.seed(42)
    GENOME_LENGTH = 4560251
    df = pd.read_csv(sample_info, sep='\t')
    SUBSAMP_COV = int(df[df['sample_id'] == sample]['cov'].values[0])

    MAX_BASES = GENOME_LENGTH * SUBSAMP_COV

    count = 0
    seqlengths = list()
    for fq in infastqs:
        for name, seq, spacer, qual in misc.fastq_reader(fq):
            seqlengths.append((hash_string(name), len(seq)))
            count += 1

    random.shuffle(seqlengths)
    total_count = 0
    keep_seqs = set()
    for s, sl in seqlengths:
        if sl < 50:
            continue
        total_count += sl
        keep_seqs.add(s)
        if total_count > MAX_BASES:
            break

    outfile = gzip.open(outfastq, 'wt')
    for fq in infastqs:
        for name, seq, spacer, qual in misc.fastq_reader(fq):
            if hash_string(name) in keep_seqs:
                print(name, seq, spacer, qual, sep='\n', file=outfile)
    outfile.close()

def extract_flanks_from_donor_reads(sample_id, infastq, sample_info, outleft, outright, min_flank_length=20):

    # Get donor sequences
    df = pd.read_csv(sample_info, sep='\t')
    df = df[df['sample_id'] == sample_id]
    LE_donor = df['le_donor'].values[0]
    RE_donor = df['re_donor'].values[0]

    # Check that donor sequences are the same length
    if len(LE_donor) != len(RE_donor):
        print("Donor lengths are not equal, exiting...")
        sys.exit(1)
    donor_end_length = len(LE_donor)

    # Get all possible edits of the donor sequence
    LE_donor_fwd_set = set()
    LE_donor_rev_set = set()
    RE_donor_fwd_set = set()
    RE_donor_rev_set = set()

    for e1 in misc.edit1(LE_donor):
        for e2 in misc.edit1(e1):
            LE_donor_fwd_set.add(e2)
            LE_donor_rev_set.add(reverse_complement(e2))

    for e1 in misc.edit1(RE_donor):
        for e2 in misc.edit1(e1):
            RE_donor_fwd_set.add(e2)
            RE_donor_rev_set.add(reverse_complement(e2))

    outleft = open(outleft, 'w')
    outright = open(outright, 'w')

    # Iterate over reads
    for read in misc.fastq_reader(infastq):
        read = list(read)
        read[1] = read[1].upper()
        for i,kmer in enumerate(misc.get_kmers(read[1], k=donor_end_length)):

            # Check if kmer is in LE_donor_fwd_set and if so, extract flanking sequence
            if kmer in LE_donor_fwd_set:
                flank_read = read[1][:i]
                flank_quals = read[3][:i]
                if len(flank_read) < min_flank_length:
                    continue
                print(read[0].split()[0]+'_leftflank', flank_read, '+', flank_quals, sep='\n', file=outleft)

            # Check if kmer is in LE_donor_rev_set and if so, extract flanking sequence
            elif kmer in LE_donor_rev_set:
                flank_read = read[1][(i+donor_end_length):]
                flank_quals = read[3][(i+donor_end_length):]
                if len(flank_read) < min_flank_length:
                    continue
                flank_read = reverse_complement(flank_read)
                flank_quals = flank_quals[::-1]
                print(read[0].split()[0]+'_leftflank', flank_read, '+', flank_quals, sep='\n', file=outleft)

            # Check if kmer is in RE_donor_fwd_set and if so, extract flanking sequence
            elif kmer in RE_donor_fwd_set:
                flank_read = read[1][i+donor_end_length:]
                flank_quals = read[3][i+donor_end_length:]
                if len(flank_read) < min_flank_length:
                    continue
                print(read[0].split()[0] + '_rightflank', flank_read, '+', flank_quals, sep='\n', file=outright)

            # Check if kmer is in RE_donor_rev_set and if so, extract flanking sequence
            elif kmer in RE_donor_rev_set:
                flank_read = read[1][:i]
                flank_quals = read[3][:i]
                if len(flank_read) < min_flank_length:
                    continue
                flank_read = reverse_complement(flank_read)
                flank_quals = flank_quals[::-1]
                print(read[0].split()[0] + '_rightflank', flank_read, '+', flank_quals, sep='\n', file=outright)

    outleft.close()
    outright.close()

def write_insertion_sites_to_tsv(insertion_sites, outfile):
    with open(outfile, 'w') as outf:
        print('genome_leftflank', 'genome_rightflank', 'plasmid_leftflank', 'plasmid_rightflank', sep='\t', file=outf)
        for j in insertion_sites:
            print(*j, sep='\t', file=outf)

def identify_insertion_sites(
        genome_leftflanks_bam, genome_rightflanks_bam, plasmid_leftflanks_bam, plasmid_rightflanks_bam,
        genome_fasta, insertion_sites_bam, insertion_sites_tsv
):

    # Load genome fasta
    genome_fasta = {rec.id:str(rec.seq) for rec in SeqIO.parse(genome_fasta, 'fasta')}

    # Load flanks aligned to genome
    # Keep only the best alignment for each read
    best_score = dict()
    for bamfile in [genome_leftflanks_bam, genome_rightflanks_bam]:
        bam = pysam.AlignmentFile(bamfile, 'rb')
        for read in bam:
            if read.is_unmapped:
                continue
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if read.mapping_quality < 20:
                continue

            ascore = float(dict(read.get_tags())['AS'])
            if read.query_name not in best_score:
                best_score[read.query_name] = (ascore, 'genome')
            elif ascore > best_score[read.query_name][0]:
                best_score[read.query_name] = (ascore, 'genome')
        bam.close()

    # Load flanks aligned to plasmid
    # Keep only the best alignment for each read
    for bamfile in [plasmid_leftflanks_bam, plasmid_rightflanks_bam]:
        bam = pysam.AlignmentFile(bamfile, 'rb')
        for read in bam:
            if read.is_unmapped:
                continue
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if read.mapping_quality < 20:
                continue

            ascore = float(dict(read.get_tags())['AS'])
            if read.query_name not in best_score:
                best_score[read.query_name] = (ascore, 'plasmid')
            elif ascore > best_score[read.query_name][0]:
                best_score[read.query_name] = (ascore, 'plasmid')
        bam.close()

    # Write insertion sites to tsv
    insertion_sites_tsv_out = open(insertion_sites_tsv, 'w')
    print("read_id", "flankside", "strand", "query_length", "contig_id", "core_start", "core_end", "genome_core",
          "insertion_core", "target_11mer", "target_14mer", "softclipped_length", sep = '\t', file = insertion_sites_tsv_out)

    # Write leftflank insertion sites to bam
    bam = pysam.AlignmentFile(genome_leftflanks_bam, 'rb')
    insertion_sites_bam_out = pysam.AlignmentFile(insertion_sites_bam, 'wb', template=bam)
    for read in bam:

        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.mapping_quality < 20:
            continue

        if read.query_name not in best_score:
            continue
        if best_score[read.query_name][1] == 'plasmid':
            continue

        strand = '+'
        if read.is_reverse:
            strand = '-'

        if strand == '+':

            right_sclen = 0
            last_pos = read.reference_end
            if sctools.is_right_softclipped_lenient(read):
                right_sclen = sctools.right_softclip_length(read)
                last_pos = sctools.right_softclipped_site_lenient(read)[1]

            if right_sclen > 10:
                continue

            core_end_pos = last_pos + right_sclen
            core_start_pos = core_end_pos - 2

            genome_core = genome_fasta[read.reference_name][core_start_pos:core_end_pos]
            insertion_core = read.query_sequence[-2:]
            target_11mer = genome_fasta[read.reference_name][core_start_pos-7:core_end_pos+2]
            target_14mer = genome_fasta[read.reference_name][core_start_pos - 7:core_end_pos + 5]

            print(read.query_name.replace("_leftflank", ""), 'LF', strand, len(read.query_sequence), read.reference_name,
                  core_start_pos, core_end_pos, genome_core, insertion_core, target_11mer, target_14mer,
                  right_sclen, sep='\t', file=insertion_sites_tsv_out)
            insertion_sites_bam_out.write(read)
        else:

            left_sclen = 0
            last_pos = read.reference_start
            if sctools.is_left_softclipped_lenient(read):
                left_sclen = sctools.left_softclip_length(read)
                last_pos = sctools.left_softclipped_site_lenient(read)[1]

            if left_sclen > 10:
                continue

            core_start_pos = last_pos
            core_end_pos = core_start_pos + 2

            genome_core = reverse_complement(genome_fasta[read.reference_name][core_start_pos:core_end_pos])
            insertion_core = reverse_complement(read.query_sequence[:2])
            target_11mer = reverse_complement(genome_fasta[read.reference_name][core_start_pos - 2:core_end_pos + 7])
            target_14mer = reverse_complement(genome_fasta[read.reference_name][core_start_pos - 5:core_end_pos + 7])

            print(read.query_name.replace("_leftflank", ""), 'LF', strand, len(read.query_sequence), read.reference_name,
                  core_start_pos, core_end_pos, genome_core, insertion_core, target_11mer, target_14mer,
                  left_sclen, sep='\t', file=insertion_sites_tsv_out)
            insertion_sites_bam_out.write(read)
    bam.close()

    # Write rightflank insertion sites to bam
    bam = pysam.AlignmentFile(genome_rightflanks_bam, 'rb')
    for read in bam:

        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.mapping_quality < 20:
            continue

        if read.query_name not in best_score:
            continue
        if best_score[read.query_name][1] == 'plasmid':
            continue

        strand = '+'
        if read.is_reverse:
            strand = '-'

        if strand == '+':

            left_sclen = 0
            last_pos = read.reference_start
            if sctools.is_left_softclipped_lenient(read):
                left_sclen = sctools.left_softclip_length(read)
                last_pos = sctools.left_softclipped_site_lenient(read)[1]

            if left_sclen > 10:
                continue

            core_start_pos = last_pos
            core_end_pos = core_start_pos + 2

            genome_core = genome_fasta[read.reference_name][core_start_pos:core_end_pos]
            insertion_core = read.query_sequence[:2]
            target_11mer = genome_fasta[read.reference_name][core_start_pos - 7:core_end_pos + 2]
            target_14mer = genome_fasta[read.reference_name][core_start_pos - 7:core_end_pos + 5]

            print(read.query_name.replace("_rightflank", ""), 'RF', strand, len(read.query_sequence), read.reference_name,
                  core_start_pos, core_end_pos, genome_core, insertion_core, target_11mer, target_14mer,
                  left_sclen, sep='\t', file=insertion_sites_tsv_out)
            insertion_sites_bam_out.write(read)
        else:

            right_sclen = 0
            last_pos = read.reference_end
            if sctools.is_right_softclipped_lenient(read):
                right_sclen = sctools.right_softclip_length(read)
                last_pos = sctools.right_softclipped_site_lenient(read)[1]

            if right_sclen > 10:
                continue

            core_end_pos = last_pos + right_sclen
            core_start_pos = core_end_pos - 2

            genome_core = reverse_complement(genome_fasta[read.reference_name][core_start_pos:core_end_pos])
            insertion_core = reverse_complement(read.query_sequence[-2:])
            target_11mer = reverse_complement(genome_fasta[read.reference_name][core_start_pos - 2:core_end_pos + 7])
            target_14mer = reverse_complement(genome_fasta[read.reference_name][core_start_pos - 5:core_end_pos + 7])

            print(read.query_name.replace("_rightflank", ""), 'RF', strand, len(read.query_sequence), read.reference_name,
                  core_start_pos, core_end_pos, genome_core, insertion_core, target_11mer, target_14mer,
                  right_sclen, sep='\t', file=insertion_sites_tsv_out)
            insertion_sites_bam_out.write(read)
    bam.close()
    insertion_sites_bam_out.close()
    insertion_sites_tsv_out.close()

    subprocess.run("samtools sort -o {} {}".format(insertion_sites_bam+'.tmp', insertion_sites_bam).split())
    os.rename(insertion_sites_bam+'.tmp', insertion_sites_bam)
    subprocess.run("samtools index {}".format(insertion_sites_bam).split())

def merge_insertion_sites(insertion_sites_tsv, merged_insertion_sites_tsv):

    insertion_sites_df = pd.read_csv(insertion_sites_tsv, sep='\t')
    if insertion_sites_df.shape[0] == 0:
        out_df = pd.DataFrame(columns=list(insertion_sites_df.columns) + ['window_start', 'window_end'])
        out_df.to_csv(merged_insertion_sites_tsv, sep='\t', index=False)
        return

    corebed = insertion_sites_df[['contig_id', 'core_start', 'core_end', 'strand']].drop_duplicates()
    corebed['Name'] = corebed['contig_id'] + ':' + corebed['core_start'].astype(str) + '-' + corebed['core_end'].astype(str)
    corebed['score'] = '.'
    corebed = corebed[['contig_id', 'core_start', 'core_end', 'Name', 'score', 'strand']]
    windowbed = pd.DataFrame.copy(corebed)
    windowbed['core_start'] = windowbed['core_start'] - 5
    windowbed['core_end'] = windowbed['core_end'] + 5

    tmpbed1 = tempfile.mktemp()
    tmpbed2 = tempfile.mktemp()
    tmpbed3 = tempfile.mktemp()

    corebed.to_csv(tmpbed1, sep='\t', header=None, index=False)
    windowbed.to_csv(tmpbed2, sep='\t', header=None, index=False)

    cmd = 'bedtools sort -i {} | bedtools merge -s -c 6 -o distinct | bedtools intersect -a stdin -b {} -wa -wb > {}'.format(
        tmpbed2, tmpbed1, tmpbed3
    )
    subprocess.run(cmd, shell=True, check=True, capture_output=True)
    merged_regions = pd.read_csv(tmpbed3, sep='\t', header=None)
    merged_regions.columns = ['contig_id', 'window_start', 'window_end', 'window_strand', 'contig_id_', 'core_start',
                              'core_end', 'core_name', 'score', 'strand']
    merged_regions = merged_regions[merged_regions.window_strand == merged_regions.strand]
    merged_regions = merged_regions[['contig_id', 'core_start', 'core_end', 'strand', 'window_start', 'window_end']]

    before = insertion_sites_df.shape[0]
    insertion_sites_df = insertion_sites_df.merge(merged_regions, on=['contig_id', 'core_start', 'core_end', 'strand'], how='inner')
    after = insertion_sites_df.shape[0]
    if before != after:
        raise Exception('Error in merging insertion site tables!')

    insertion_sites_df.to_csv(merged_insertion_sites_tsv, sep='\t', index=False)

def remove_plasmid_alignments(genome_bam, plasmid_bam, outbam):

    plasmid_reads = set()
    for i,read in enumerate(pysam.AlignmentFile(plasmid_bam, 'rb')):
        if read.is_unmapped:
            continue
        plasmid_reads.add(read.query_name)

    genome_bam = pysam.AlignmentFile(genome_bam, 'rb')
    outbam = pysam.AlignmentFile(outbam, 'wb', template=genome_bam)
    for read in genome_bam:
        if read.is_unmapped:
            continue
        if read.query_name in plasmid_reads:
            continue
        outbam.write(read)
    outbam.close()

    cmd = ['samtools', 'index', outbam.filename]
    subprocess.run(cmd, check=True)

def extract_indels_from_read(read, minsize=50):
    # Extract all indels from read using pysam

    qry_pos = 0
    ref_pos = read.reference_start
    ref_read_pos = 0
    indels = []
    for op, length in read.cigartuples:
        if op == 0:
            qry_pos += length
            ref_pos += length
            ref_read_pos += length
        elif op == 1:
            if length >= minsize:
                indels.append(('I', qry_pos, ref_pos, length, read.query_sequence[qry_pos:qry_pos+length]))
            qry_pos += length
        elif op == 2:
            if length >= minsize:
                indels.append(('D', qry_pos, ref_pos, length, read.get_reference_sequence()[ref_read_pos:ref_read_pos+length]))
            ref_pos += length
            ref_read_pos += length
        elif op == 4:
            qry_pos += length
        elif op == 5:
            pass
        else:
            print("Unexpected CIGAR operation")
            print(op, length)
            sys.exit(1)

    return indels

def extract_indels(bam, outfile):
    outfile = open(outfile, 'w')
    print("read_name", "ref_name", "read_ref_start", "read_ref_end", "read_strand", "indel_type", "qry_pos", "ref_pos", "length", "seq", sep='\t', file=outfile)
    for read in pysam.AlignmentFile(bam, 'rb'):
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        strand = '-' if read.is_reverse else '+'

        if "D" not in read.cigarstring and "I" not in read.cigarstring:
            continue

        has_longdel = False
        has_longins = False
        if "D" in read.cigarstring:
            for dellen in misc.get_deletion_lengths(read.cigarstring):
                if dellen >= 50:
                    has_longdel = True
                    break
        if "I" in read.cigarstring:
            for inslen in misc.get_insertion_lengths(read.cigarstring):
                if inslen >= 50:
                    has_longins = True
                    break

        if not has_longdel and not has_longins:
            continue

        indels = extract_indels_from_read(read, 50)
        for indel in indels:
            print(read.query_name, read.reference_name, read.reference_start, read.reference_end, strand, *indel, sep='\t', file=outfile)
    outfile.close()

def cluster_indels(indel_files, outtsv, threads):

    # Merge all indel files into one temporary file
    outtmp = open(outtsv+'.tmp', 'w')
    print('sample_id', 'biorep', 'read_name', 'ref_name', 'read_ref_start', 'read_ref_end',
          'read_strand', 'indel_type', 'qry_pos', 'ref_pos', 'length', 'seq', sep='\t', file=outtmp)
    for f in indel_files:
        sample_id, biorep = os.path.basename(f).split('.')[:2]
        with open(f) as infile:
            infile.readline()
            for line in infile:
                print(sample_id, biorep, *line.strip().split('\t'), sep='\t', file=outtmp)
    outtmp.close()

    # Write indels to bed and fasta
    tmpfasta = open(outtmp.name + '.fasta', 'w')
    tmpbed = open(outtmp.name + '.bed', 'w')
    with open(outtmp.name) as infile:
        header = infile.readline().strip().split('\t')
        for i, line in enumerate(infile):
            line = line.strip().split('\t')
            line = dict(zip(header, line))
            print('>' + str(i), line['seq'], sep='\n', file=tmpfasta)
            indel_name = line['read_name']+':'+line['ref_name']+':'+line['read_strand']+'_'+line['indel_type']+'_'+line['qry_pos']+'_'+line['ref_pos']
            print(line['ref_name'], int(line['ref_pos'])-20, int(line['ref_pos'])+20, indel_name, sep='\t', file=tmpbed)
    tmpfasta.close()
    tmpbed.close()

    outpref = outtsv.replace(".tsv", "") + '.mmseqs'
    # Cluster indels with mmseqs
    cmd = '{mmseqs} easy-cluster --threads {threads} --cov-mode 0 --min-seq-id 0.9 -c 0.9 {infasta} {outprefix} {tmpdir}'.format(
        mmseqs=MMSEQS2_PATH, threads=threads, infasta=tmpfasta.name, outprefix=outpref, tmpdir='tmp'
    )
    subprocess.run(cmd, shell=True, check=True)
    shutil.rmtree('tmp')
    cluster_assignment = dict()
    with open(outpref + '_cluster.tsv') as infile:
        for line in infile:
            line = line.strip().split('\t')
            cluster_assignment[int(line[1])] = int(line[0])
    os.remove(tmpfasta.name)

    # Merge indels in the same cluster
    tmpbedmerged = tmpbed.name+'.merged'
    cmd = 'bedtools sort -i {inbed} | bedtools merge -i stdin -c 4 -o distinct > {out}'.format(inbed=tmpbed.name, out=tmpbedmerged)
    subprocess.run(cmd, shell=True, check=True)
    name2merged = dict()
    with open(tmpbedmerged) as infile:
        for line in infile:
            line = line.strip().split('\t')
            merged = line[0]+':'+line[1]+'-'+line[2]
            for name in line[3].split(','):
                name2merged[name] = merged
    os.remove(tmpbed.name)
    os.remove(tmpbedmerged)

    # Write output tsv
    outtsv = open(outtsv, 'w')
    with open(outtmp.name) as infile:
        header = infile.readline().strip().split('\t')
        print(*header, 'seq_cluster', 'merged_position', sep='\t', file=outtsv)
        for i,line in enumerate(infile):
            line_orig = line.strip().split('\t')
            line = dict(zip(header, line_orig))
            indel_name = line['read_name'] + ':' + line['ref_name'] + ':' + line['read_strand'] + '_' + line['indel_type'] + '_' + line['qry_pos'] + '_' + line['ref_pos']
            print(*line.values(), cluster_assignment[i], name2merged[indel_name], sep='\t', file=outtsv)
    outtsv.close()
    os.remove(outtmp.name)

def fgsv_filter_size(intxt, genome_fasta, outtxt):
    # Filter out breakpoints that are below 50 bp in length
    genome_len = len([rec for rec in SeqIO.parse(genome_fasta, 'fasta')][0].seq)
    df = pd.read_csv(intxt, sep='\t')
    filtered_df = []
    for i, row in df.iterrows():
        left_pos = row['left_pos']
        right_pos = row['right_pos']

        if left_pos == 1 or right_pos == genome_len:
            continue

        noncircular_dist = right_pos-(left_pos-1)
        circular_dist = genome_len - noncircular_dist
        if noncircular_dist < 50 or circular_dist < 50:
            continue
        filtered_df.append(row)

    filtered_df = pd.DataFrame(filtered_df).reset_index(drop=True)
    filtered_df.to_csv(outtxt, sep='\t', index=False)

def merge_breakpoints(intsvs, outtsv):

    left_breakpoint_tempbedfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    right_breakpoint_tempbedfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for tsv in intsvs:
        with open(tsv) as infile:
            header = infile.readline().strip().split('\t')
            for line in infile:
                line = line.strip().split('\t')
                line = dict(zip(header, line))
                left_name = line['left_contig'] + ':' + line['left_min_pos'] + ':' + line['left_strand']
                right_name = line['right_contig'] + ':' + line['right_min_pos'] + ':' + line['right_strand']

                left_start_pos = int(line['left_min_pos']) - 20
                left_end_pos = int(line['left_min_pos']) + 20
                right_start_pos = int(line['right_min_pos']) - 20
                right_end_pos = int(line['right_min_pos']) + 20

                if left_start_pos < 0:
                    left_start_pos = 0
                if right_start_pos < 0:
                    right_start_pos = 0

                print(line['left_contig'], left_start_pos, left_end_pos, left_name, '.', line['left_strand'],
                      sep='\t', file=left_breakpoint_tempbedfile)
                print(line['right_contig'], right_start_pos, right_end_pos, right_name, '.', line['right_strand'],
                        sep='\t', file=right_breakpoint_tempbedfile)
    left_breakpoint_tempbedfile.close()
    right_breakpoint_tempbedfile.close()

    cmd = "bedtools sort -i {} | bedtools merge -s -c 4,6 -o distinct > {}".format(left_breakpoint_tempbedfile.name, left_breakpoint_tempbedfile.name + '.merged')
    subprocess.run(cmd, shell=True, check=True)
    cmd = "bedtools sort -i {} | bedtools merge -s -c 4,6 -o distinct > {}".format(right_breakpoint_tempbedfile.name, right_breakpoint_tempbedfile.name + '.merged')
    subprocess.run(cmd, shell=True, check=True)

    left_name2merged = dict()
    with open(left_breakpoint_tempbedfile.name+'.merged') as infile:
        for line in infile:
            line = line.strip().split('\t')
            merged = line[0] + ':' + line[1] + '-' + line[2] + ':' + line[4]
            for name in line[3].split(','):
                left_name2merged[name] = merged

    right_name2merged = dict()
    with open(right_breakpoint_tempbedfile.name+'.merged') as infile:
        for line in infile:
            line = line.strip().split('\t')
            merged = line[0] + ':' + line[1] + '-' + line[2] + ':' + line[4]
            for name in line[3].split(','):
                right_name2merged[name] = merged

    os.remove(left_breakpoint_tempbedfile.name)
    os.remove(left_breakpoint_tempbedfile.name + '.merged')
    os.remove(right_breakpoint_tempbedfile.name)
    os.remove(right_breakpoint_tempbedfile.name + '.merged')

    wrote_header = False
    with open(outtsv, 'w') as outfile:
        for tsv in intsvs:
            fname = os.path.basename(tsv).split(".")
            sample_id = fname[0]
            biorep = fname[1]
            with open(tsv) as infile:
                header = infile.readline().strip().split('\t')
                if not wrote_header:
                    print("sample_id", "biorep", '\t'.join(header[:-2]), "left_merged_breakpoint", "right_merged_breakpoint", sep='\t', file=outfile)
                    wrote_header = True
                for line in infile:
                    line_orig = line.strip().split('\t')
                    line = dict(zip(header, line_orig))
                    left_name = line['left_contig'] + ':' + line['left_min_pos'] + ':' + line['left_strand']
                    right_name = line['right_contig'] + ':' + line['right_min_pos'] + ':' + line['right_strand']
                    print(sample_id, biorep, '\t'.join(line_orig), left_name2merged[left_name], right_name2merged[right_name], sep='\t', file=outfile)

def overlapping_breakpoints_indels(indels_tsv, breakpoints_tsv, outtsv):

    indels = pd.read_csv(indels_tsv, sep='\t')
    breakpoints = pd.read_csv(breakpoints_tsv, sep='\t')

    uniq_regions = set()
    for merged_position, indel_type, medseqlen in zip(indels.merged_position, indels.indel_type, indels.medseqlen):
        chrom, pos = merged_position.split(':')
        start, end = pos.split('-')
        if indel_type == 'I':
            uniq_regions.add((chrom, start, end, merged_position+':'+indel_type))
        elif indel_type == 'D':
            uniq_regions.add((chrom, start, end, merged_position+':'+indel_type+"_left"))
            uniq_regions.add((chrom, int(start)+int(medseqlen), int(end)+int(medseqlen), merged_position+':'+indel_type+"_right"))

    indel_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for chrom, start, end, merged_position in uniq_regions:
        print(chrom, start, end, merged_position, sep='\t', file=indel_tmpbed)
    indel_tmpbed.close()

    uniq_regions = set()
    for left_merged_breakpoint, right_merged_breakpoint in zip(breakpoints.left_merged_breakpoint, breakpoints.right_merged_breakpoint):
        chrom, pos, strand = left_merged_breakpoint.split(':')
        start, end = pos.split('-')
        uniq_regions.add((chrom, start, end, left_merged_breakpoint))
        chrom, pos, strand = right_merged_breakpoint.split(':')
        start, end = pos.split('-')
        uniq_regions.add((chrom, start, end, right_merged_breakpoint))

    breakpoints_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for chrom, start, end, merged_position in uniq_regions:
        print(chrom, start, end, merged_position, sep='\t', file=breakpoints_tmpbed)
    breakpoints_tmpbed.close()

    intersect_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    intersect_tmpbed.close()
    cmd = "bedtools intersect -a {} -b {} -wa -wb > {}".format(indel_tmpbed.name, breakpoints_tmpbed.name, intersect_tmpbed.name)
    subprocess.run(cmd, shell=True, check=True)
    outtsv = open(outtsv, 'w')
    print("indel", "breakpoint", sep='\t', file=outtsv)
    with open(intersect_tmpbed.name) as infile:
        for line in infile:
            line = line.strip().split('\t')
            indel = line[3]
            breakpoint = line[7]
            print(indel, breakpoint, sep='\t', file=outtsv)
    outtsv.close()

    os.remove(indel_tmpbed.name)
    os.remove(breakpoints_tmpbed.name)
    os.remove(intersect_tmpbed.name)

def overlapping_insertion_sites(indels_tsv, breakpoints_tsv, programmed_insertions_tsv, wt_insertions_tsv, outtsv):

    indels = pd.read_csv(indels_tsv, sep='\t')
    breakpoints = pd.read_csv(breakpoints_tsv, sep='\t')
    programmed_insertions = pd.read_csv(programmed_insertions_tsv, sep='\t')
    wt_insertions = pd.read_csv(wt_insertions_tsv, sep='\t')

    uniq_regions = set()
    for merged_position, indel_type, medseqlen in zip(indels.merged_position, indels.indel_type, indels.medseqlen):
        chrom, pos = merged_position.split(':')
        start, end = pos.split('-')
        if indel_type == 'I':
            uniq_regions.add((chrom, start, end, merged_position + ':' + indel_type))
        elif indel_type == 'D':
            uniq_regions.add((chrom, start, end, merged_position + ':' + indel_type + "_left"))
            uniq_regions.add((chrom, int(start) + int(medseqlen), int(end) + int(medseqlen),
                              merged_position + ':' + indel_type + "_right"))

    indel_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for chrom, start, end, merged_position in uniq_regions:
        print(chrom, start, end, merged_position, sep='\t', file=indel_tmpbed)
    indel_tmpbed.close()

    uniq_regions = set()
    for left_merged_breakpoint, right_merged_breakpoint in zip(breakpoints.left_merged_breakpoint,
                                                               breakpoints.right_merged_breakpoint):
        chrom, pos, strand = left_merged_breakpoint.split(':')
        start, end = pos.split('-')
        uniq_regions.add((chrom, start, end, left_merged_breakpoint))
        chrom, pos, strand = right_merged_breakpoint.split(':')
        start, end = pos.split('-')
        uniq_regions.add((chrom, start, end, right_merged_breakpoint))

    breakpoints_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for chrom, start, end, merged_position in uniq_regions:
        print(chrom, start, end, merged_position, sep='\t', file=breakpoints_tmpbed)
    breakpoints_tmpbed.close()

    insertion_sites_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for bridgerna_id, rtg_extension, contig_id, core_start, core_end, target_11mer, exp_target_11mer, offtarget_category in zip(
            programmed_insertions.bridgerna_id, programmed_insertions.rtg_extension,
            programmed_insertions.contig_id, programmed_insertions.core_start, programmed_insertions.core_end,
            programmed_insertions.target_11mer, programmed_insertions.exp_target_11mer,
            programmed_insertions.offtarget_category):
        name = contig_id+":"+str(core_start)+":"+bridgerna_id.replace(" ", "-")+":"+rtg_extension.replace(" ", "-")+\
               ":"+target_11mer+':'+exp_target_11mer+':'+offtarget_category.replace(" ", "-")
        print(contig_id, core_start, core_end, name, sep='\t', file=insertion_sites_tmpbed)
    for contig_id, core_start, core_end, target_11mer, exp_target_11mer in zip(
            wt_insertions.contig_id, wt_insertions.core_start, wt_insertions.core_end,
            wt_insertions.target_11mer, wt_insertions.exp_target_11mer):
        name = contig_id+":"+str(core_start)+":WT:"+target_11mer+':'+exp_target_11mer
        print(contig_id, core_start, core_end, name, sep='\t', file=insertion_sites_tmpbed)
    insertion_sites_tmpbed.close()

    intersect_tmpbed1 = tempfile.NamedTemporaryFile(mode='w', delete=False)
    intersect_tmpbed1.close()
    cmd = "bedtools intersect -a {} -b {} -wa -wb > {}".format(indel_tmpbed.name, insertion_sites_tmpbed.name, intersect_tmpbed1.name)
    subprocess.run(cmd, shell=True, check=True)

    intersect_tmpbed2 = tempfile.NamedTemporaryFile(mode='w', delete=False)
    intersect_tmpbed2.close()
    cmd = "bedtools intersect -a {} -b {} -wa -wb > {}".format(breakpoints_tmpbed.name, insertion_sites_tmpbed.name, intersect_tmpbed2.name)
    subprocess.run(cmd, shell=True, check=True)

    with open(outtsv, 'w') as outtsv:
        print("strucvar_type", "strucvar_position", "insertion_name", sep='\t', file=outtsv)
        with open(intersect_tmpbed1.name) as infile:
            for line in infile:
                line = line.strip().split('\t')
                insertion = line[7]
                indel = line[3]
                print("indel", indel, insertion, sep='\t', file=outtsv)
        with open(intersect_tmpbed2.name) as infile:
            for line in infile:
                line = line.strip().split('\t')
                insertion = line[7]
                breakpoint = line[3]
                print("breakpoint", breakpoint, insertion, sep='\t', file=outtsv)

    os.remove(indel_tmpbed.name)
    os.remove(breakpoints_tmpbed.name)
    os.remove(insertion_sites_tmpbed.name)
    os.remove(intersect_tmpbed1.name)
    os.remove(intersect_tmpbed2.name)

def indel_site_coverage(indels, inbams, threads, outtsv):

    indels = pd.read_csv(indels, sep='\t')

    uniq_regions = set()
    for sample_id, indel_position, in zip(indels.sample_id, indels.indel_position):
        chrom, pos = indel_position.split(':')
        start = int(pos)-1
        end = start + 1
        uniq_regions.add((chrom, start, end, indel_position))

    indel_tmpbed = tempfile.NamedTemporaryFile(mode='w', delete=False)
    for region in uniq_regions:
        print(*region, sep='\t', file=indel_tmpbed)
    indel_tmpbed.close()

    with open(outtsv, 'w') as outtsv:
        print("sample_id", "biorep", "indel_position", "coverage", sep='\t', file=outtsv)
        for inbam in inbams:
            print("Processing", inbam)
            sample_id, biorep = os.path.basename(inbam).split('.')[:2]
            tmp_outtsv = tempfile.NamedTemporaryFile(mode='w', delete=False)
            tmp_outtsv.close()
            cmd = "samtools depth --threads {} -b {} {} > {}".format(threads, indel_tmpbed.name, inbam, tmp_outtsv.name)
            subprocess.run(cmd, shell=True, check=True)
            with open(tmp_outtsv.name) as infile:
                for line in infile:
                    line = line.strip().split('\t')
                    print(sample_id, biorep, line[0]+':'+line[1], line[2], sep='\t', file=outtsv)
            os.remove(tmp_outtsv.name)

    os.remove(indel_tmpbed.name)
















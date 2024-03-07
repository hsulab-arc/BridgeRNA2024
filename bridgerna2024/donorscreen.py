import subprocess
from bridgerna2024 import *
from bridgerna2024 import misc
import pysam, os
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import reverse_complement

def _run_donorscreen_snakemake(workdir, threads):
    """Run the RNA-seq snakemake pipeline."""
    subprocess.run(['snakemake', '--snakefile', DONORSCREEN_SNAKEFILE_PATH, '--configfile',
                    DONORSCREEN_SNAKEMAKE_CONFIG_PATH, '-j', str(threads),
                    '--config',  'workdir={}'.format(workdir)], check=True)

def umi_align_extract_variable_regions(bampath, amplicon_fasta, amplicon_bed, outdir):

    # Extract sample name from file path
    sample = os.path.basename(bampath).split('.')[0]

    # Get amplicon sequence
    amplicon = None
    for rec in SeqIO.parse(amplicon_fasta, 'fasta'):
        amplicon = str(rec.seq)
    amp_length = len(amplicon)

    # Specify the different variable regions of the amplicon
    region_names = ['donor_m10box', 'm35box', 'ldg', 'rdg', 'ltg_core', 'rtg_core', 'umi']

    # Extract the variable regions from the amplicon bed file
    position2region = dict()
    region2position = dict()
    with open(amplicon_bed) as infile:
        for line in infile:
            line = line.strip().split('\t')
            if line[3] not in region_names:
                continue

            start, end = int(line[1]), int(line[2])

            for i in range(start, end):
                position2region[i] = line[3]
            region2position[line[3]] = list(range(start, end))

    tsv_out = open(os.path.join(outdir, '{}.variable_regions.tsv'.format(sample)), 'w')
    print("read_hash_id", 'umi', 'donor_m10box', 'm35box', 'ldg', 'rdg', 'ltg_core', 'rtg_core', sep='\t', file=tsv_out)

    bam = pysam.AlignmentFile(bampath, 'rb')
    for reads1, reads2 in misc.retrieve_read_pairs_from_name_sorted_bam(bam):
        read1, read2 = None, None
        for read in reads1:
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            # Skip mate unmapped reads
            if read.mate_is_unmapped:
                continue
            # Skip secondary alignments
            if read.is_secondary:
                continue
            # Skip supplementary alignments
            if read.is_supplementary:
                continue
            # Skip reads with wrong insert size
            readtlen = abs(read.tlen)
            if readtlen < amp_length - 5 or readtlen > amp_length + 5:
                continue
            read1 = read

        for read in reads2:
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            # Skip mate unmapped reads
            if read.mate_is_unmapped:
                continue
            # Skip secondary alignments
            if read.is_secondary:
                continue
            # Skip supplementary alignments
            if read.is_supplementary:
                continue
            # Skip reads with wrong insert size
            readtlen = abs(read.tlen)
            if readtlen < amp_length - 5 or readtlen > amp_length + 5:
                continue
            read2 = read

        if read1 is None or read2 is None:
            continue

        read1_query_seq = read1.query_sequence
        read1_query_quals = read1.query_qualities
        read2_query_seq = read2.query_sequence
        read2_query_quals = read2.query_qualities

        read1_info = {'pos2nuc': dict(), 'pos2qual': dict()}
        read2_info = {'pos2nuc': dict(), 'pos2qual': dict()}
        read1_pairs = read1.get_aligned_pairs()
        read2_pairs = read2.get_aligned_pairs()

        for query_pos, ref_pos in read1_pairs:
            if query_pos is None:
                continue
            if ref_pos not in position2region:
                continue
            try:
                read1_info['pos2nuc'][ref_pos] = read1_query_seq[query_pos]
                read1_info['pos2qual'][ref_pos] = read1_query_quals[query_pos]
            except Exception as e:
                print(query_pos)
                raise e

        for query_pos, ref_pos in read2_pairs:
            if query_pos is None:
                continue
            if ref_pos not in position2region:
                continue
            try:
                read2_info['pos2nuc'][ref_pos] = read2_query_seq[query_pos]
                read2_info['pos2qual'][ref_pos] = read2_query_quals[query_pos]
            except Exception as e:
                print(query_pos)
                raise e

        readhash = misc.hash_string(read1.query_name)

        region_consenses = dict()
        for region in region_names:

            region_consensus = ''
            for rpos in region2position[region]:

                if rpos not in read1_info['pos2nuc'] and rpos not in read2_info['pos2nuc']:
                    region_consensus += 'N'
                else:
                    r1_nuc, r1_qual = "N", 0
                    r2_nuc, r2_qual = "N", 0
                    if rpos in read1_info['pos2nuc']:
                        r1_nuc, r1_qual = read1_info['pos2nuc'][rpos], read1_info['pos2qual'][rpos]
                    if rpos in read2_info['pos2nuc']:
                        r2_nuc, r2_qual = read2_info['pos2nuc'][rpos], read2_info['pos2qual'][rpos]

                    if r1_qual < 20 and r2_qual < 20:
                        region_consensus += 'N'
                    elif r1_qual >= r2_qual:
                        region_consensus += r1_nuc
                    else:
                        region_consensus += r2_nuc

            region_consenses[region] = region_consensus

        print(readhash, region_consenses['umi'], region_consenses['donor_m10box'], region_consenses['m35box'],
              region_consenses['ldg'], region_consenses['rdg'], region_consenses['ltg_core'],
              region_consenses['rtg_core'], sep='\t', file=tsv_out)

    bam.close()
    tsv_out.close()

def map_oligos(input_umis, oligo_info, outfile):

    # Hard coded, comes from the oligo design
    OLIGO_REGIONS = {}
    OLIGO_REGIONS['M35BOX'] = (55, 59)
    OLIGO_REGIONS['DONOR_M10BOX'] = (67, 84)
    OLIGO_REGIONS['LDG'] = (212, 221)
    OLIGO_REGIONS['RDG'] = (251, 255)
    OLIGO_REGIONS['LTG_CORE'] = (146, 148)
    OLIGO_REGIONS['RTG_CORE'] = (168, 170)

    # Load the oligo design information
    oligo_map = defaultdict(set)
    with open(oligo_info) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            line = dict(zip(header, line))

            oligo = line['oligo']
            donor_m10box = oligo[OLIGO_REGIONS['DONOR_M10BOX'][0]:OLIGO_REGIONS['DONOR_M10BOX'][1]]
            m35box = oligo[OLIGO_REGIONS['M35BOX'][0]:OLIGO_REGIONS['M35BOX'][1]]
            ldg = oligo[OLIGO_REGIONS['LDG'][0]:OLIGO_REGIONS['LDG'][1]]
            rdg = oligo[OLIGO_REGIONS['RDG'][0]:OLIGO_REGIONS['RDG'][1]]
            ltg_core = oligo[OLIGO_REGIONS['LTG_CORE'][0]:OLIGO_REGIONS['LTG_CORE'][1]]
            rtg_core = oligo[OLIGO_REGIONS['RTG_CORE'][0]:OLIGO_REGIONS['RTG_CORE'][1]]

            oligo_map[(donor_m10box, m35box, ldg, rdg, ltg_core, rtg_core)].add(line['oligo_id'])

    o = open(outfile, 'w')
    print("umi", "oligo_id", sep='\t', file=o)
    with open(input_umis) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            line = dict(zip(header, line))
            if 'N' in line['umi']+line['donor_m10box']+line['m35box']+line['ldg']+line['rdg']+line['ltg_core']+line['rtg_core']:
                continue
            key = (line['donor_m10box'], line['m35box'], line['ldg'], line['rdg'], line['ltg_core'], line['rtg_core'])
            if key not in oligo_map:
                continue
            oligos = oligo_map[key]
            for oligo in oligos:
                print(line['umi'], oligo, sep='\t', file=o)
    o.close()

def final_umi_mapping(umi_files, oligos_info, final):

    oligo_donor_loops = dict()
    with open(oligos_info) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            line = dict(zip(header, line))
            ldg, rdg = line['oligo'][212:220], line['oligo'][251:255]
            oligo_donor_loops[line['oligo_id']] = (ldg, rdg)

    A_HOMOPOL, C_HOMOPOL, G_HOMOPOL, T_HOMOPOL = 'A'*6, 'C'*6, 'G'*6, 'T'*6

    umi_counts = defaultdict(int)
    for f in umi_files:
        with open(f) as infile:
            header = infile.readline().strip().split('\t')
            for line in infile:
                line = line.strip().split('\t')
                line = dict(zip(header, line))
                if A_HOMOPOL in line['umi'] or C_HOMOPOL in line['umi'] or G_HOMOPOL in line['umi'] or T_HOMOPOL in line['umi']:
                    continue
                ldg, rdg = oligo_donor_loops[line['oligo_id']]
                umi, ldg, rdg = reverse_complement(line['umi']), reverse_complement(ldg), reverse_complement(rdg)
                umi = umi + ldg + rdg
                umi_counts[(line['oligo_id'], umi)] += 1

    # Remove umis that map to multiple oligos
    umioligo_counts = defaultdict(int)
    for k in umi_counts:
        umioligo_counts[k[1]] += 1
    keep_umis = set()

    for k in umioligo_counts:
        if umioligo_counts[k] == 1:
            keep_umis.add(k)

    # Write out the final umi mapping
    o = open(final, 'w')
    print("umi", "oligo_id", 'umi_count', sep='\t', file=o)
    for k in umi_counts:
        if k[1] not in keep_umis:
            continue
        print(k[1], k[0], umi_counts[k], sep='\t', file=o)
    o.close()

def screen_extract_umi(bampath, amplicon_fasta, amplicon_bed, outtsv):

    library = os.path.basename(bampath).split('.')[0]

    # Load the amplicon information
    amp_length = None
    for rec in SeqIO.parse(amplicon_fasta, 'fasta'):
        amp_length = len(rec.seq)

    # Load the amplicon bed information
    pos2region = dict()
    barcode_range = []
    with open(amplicon_bed) as infile:
        for line in infile:
            line = line.strip().split('\t')
            rng = list(range(int(line[1]), int(line[2])))
            for p in rng:
                pos2region[p] = line[3]
            barcode_range += rng

    tsv_out = open(outtsv, 'w')
    print("read_id", 'r1_umi', 'r2_umi', 'merged_umi', sep='\t', file=tsv_out)

    bam = pysam.AlignmentFile(bampath, 'rb')
    for reads1, reads2 in misc.retrieve_read_pairs_from_name_sorted_bam(bam):
        read1, read2 = None, None
        for read in reads1:
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            # Skip mate unmapped reads
            if read.mate_is_unmapped:
                continue
            # Skip secondary alignments
            if read.is_secondary:
                continue
            # Skip supplementary alignments
            if read.is_supplementary:
                continue
            # Skip reads with wrong insert size
            readtlen = abs(read.tlen)
            if readtlen < amp_length - 5 or readtlen > amp_length + 5:
                continue
            read1 = read

        for read in reads2:
            # Skip unmapped reads
            if read.is_unmapped:
                continue
            # Skip mate unmapped reads
            if read.mate_is_unmapped:
                continue
            # Skip secondary alignments
            if read.is_secondary:
                continue
            # Skip supplementary alignments
            if read.is_supplementary:
                continue
            # Skip reads with wrong insert size
            readtlen = abs(read.tlen)
            if readtlen < amp_length - 5 or readtlen > amp_length + 5:
                continue
            read2 = read

        if read1 is None or read2 is None:
            continue

        read1_query_seq = read1.query_sequence
        read1_query_quals = read1.query_qualities
        read2_query_seq = read2.query_sequence
        read2_query_quals = read2.query_qualities

        read1_info = {'pos2nuc': dict(), 'pos2qual': dict()}
        read2_info = {'pos2nuc': dict(), 'pos2qual': dict()}
        read1_pairs = read1.get_aligned_pairs()
        read2_pairs = read2.get_aligned_pairs()

        for query_pos, ref_pos in read1_pairs:
            if query_pos is None:
                continue
            if ref_pos not in barcode_range:
                continue
            try:
                read1_info['pos2nuc'][ref_pos] = read1_query_seq[query_pos]
                read1_info['pos2qual'][ref_pos] = read1_query_quals[query_pos]
            except Exception as e:
                print(query_pos)
                raise e

        for query_pos, ref_pos in read2_pairs:
            if query_pos is None:
                continue
            if ref_pos not in barcode_range:
                continue
            try:
                read2_info['pos2nuc'][ref_pos] = read2_query_seq[query_pos]
                read2_info['pos2qual'][ref_pos] = read2_query_quals[query_pos]
            except Exception as e:
                print(query_pos)
                raise e

        r1 = read1_info
        r2 = read2_info

        r1_regions = {r:'' for r in list(set(pos2region.values()))}
        r2_regions = {r: '' for r in list(set(pos2region.values()))}
        merged_regions = {r: '' for r in list(set(pos2region.values()))}

        for pos in barcode_range:
            if pos not in r1['pos2nuc']:
                r1_regions[pos2region[pos]] += 'N'

            else:
                r1_regions[pos2region[pos]] += r1['pos2nuc'][pos]

            if pos not in r2['pos2nuc']:
                r2_regions[pos2region[pos]] += 'N'
            else:
                r2_regions[pos2region[pos]] += r2['pos2nuc'][pos]

            if pos not in r1['pos2qual'] and pos not in r2['pos2qual']:
                merged_regions[pos2region[pos]] += 'N'
            elif pos in r1['pos2qual'] and pos not in r2['pos2qual']:
                merged_regions[pos2region[pos]] += r1['pos2nuc'][pos]
            elif pos in r2['pos2qual'] and pos not in r1['pos2qual']:
                merged_regions[pos2region[pos]] += r2['pos2nuc'][pos]
            elif r1['pos2qual'][pos] >= r2['pos2qual'][pos]:
                merged_regions[pos2region[pos]] += r1['pos2nuc'][pos]
            else:
                merged_regions[pos2region[pos]] += r2['pos2nuc'][pos]

        r1_umi = r1_regions['umi'] + r1_regions['ldg'] + r1_regions['rdg']
        r2_umi = r2_regions['umi'] + r2_regions['ldg'] + r2_regions['rdg']
        merged_umi = merged_regions['umi'] + merged_regions['ldg'] + merged_regions['rdg']

        print(read1.query_name, r1_umi, r2_umi, merged_umi, sep='\t', file=tsv_out)

    bam.close()
    tsv_out.close()

def screen_umi_counts(screen_umis, umi_map, outcounts):

    # Load umi map
    umi2oligos = dict()
    with open(umi_map) as infile:
        header = infile.readline().strip().split('\t')
        for line in infile:
            line = line.strip().split('\t')
            line = dict(zip(header, line))
            umi2oligos[line['umi']] = line['oligo_id']

    # Load umi counts
    umi_counts = defaultdict(lambda: defaultdict(int))
    for f in screen_umis:
        sample = os.path.basename(f).split('.')[0]
        with open(f) as infile:

            header = infile.readline().strip().split('\t')
            for line in infile:
                line = line.strip().split('\t')
                line = dict(zip(header, line))
                if line['merged_umi'] not in umi2oligos:
                    continue
                umi_counts[sample][umi2oligos[line['merged_umi']]] += 1

    # Write out counts
    with open(outcounts, 'w') as out:
        print("sample_id", "oligo_id", "oligo_count", sep='\t', file=out)
        for samp in umi_counts:
            for oligo in umi_counts[samp]:
                print(samp, oligo, umi_counts[samp][oligo], sep='\t', file=out)
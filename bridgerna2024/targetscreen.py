import subprocess
from bridgerna2024 import *
from bridgerna2024 import misc
import pysam, re, os
from collections import defaultdict
from Bio import SeqIO

def _run_targetscreen_snakemake(workdir, threads):
    """Run the RNA-seq snakemake pipeline."""
    subprocess.run(['snakemake', '--snakefile', TARGETSCREEN_SNAKEFILE_PATH, '--configfile',
                    TARGETSCREEN_SNAKEMAKE_CONFIG_PATH, '-j', str(threads),
                    '--config',  'workdir={}'.format(workdir)], check=True)

def extract_barcodes(bampath, amplicon_fasta, outtsv):

    # Get sample name
    sample = os.path.basename(bampath).split('.')[0]

    # Get amplicon sequence
    amplicon = None
    for rec in SeqIO.parse(amplicon_fasta, 'fasta'):
        amplicon = str(rec.seq)
    amp_length = len(amplicon)

    # Get barcode range
    barcode_start, barcode_end = re.search(r"N{12}", amplicon).span()
    barcode_range = list(range(barcode_start, barcode_end))

    # Get barcode info
    read_info = defaultdict(lambda : {"R1":{'pos2nuc':{}, 'pos2qual': {}},
                                      "R2":{'pos2nuc':{}, 'pos2qual':{}}})

    # Write barcode info to tsv
    tsv_out = open(outtsv, 'w')
    print("read_id", "r1_barcode", "r2_barcode", "merged_barcode", sep='\t', file=tsv_out)
    bam = pysam.AlignmentFile(bampath, "rb")
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
            # Skip reads with low mapping quality
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
            # Skip reads with incorrect insert size
            readtlen = abs(read.tlen)
            if readtlen < amp_length - 5 or readtlen > amp_length + 5:
                continue
            read2 = read

        if read1 is None or read2 is None:
            continue

        # Get barcode sequence and quality
        read1_query_seq = read1.query_sequence
        read1_query_quals = read1.query_qualities
        read2_query_seq = read2.query_sequence
        read2_query_quals = read2.query_qualities

        read1_info = {'pos2nuc': dict(), 'pos2qual':dict()}
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

        # Get barcode sequences
        r1 = read1_info
        r2 = read2_info

        r1_barcode = ""
        r2_barcode = ""
        merged_barcode = ""

        for pos in barcode_range:
            if pos not in r1['pos2nuc']:
                r1_barcode += 'N'
            else:
                r1_barcode += r1['pos2nuc'][pos]

            if pos not in r2['pos2nuc']:
                r2_barcode += 'N'
            else:
                r2_barcode += r2['pos2nuc'][pos]

            if pos not in r1['pos2qual'] and pos not in r2['pos2qual']:
                merged_barcode += 'N'
            elif pos in r1['pos2qual'] and pos not in r2['pos2qual']:
                merged_barcode += r1['pos2nuc'][pos]
            elif pos in r2['pos2qual'] and pos not in r1['pos2qual']:
                merged_barcode += r2['pos2nuc'][pos]
            elif r1['pos2qual'][pos] >= r2['pos2qual'][pos]:
                merged_barcode += r1['pos2nuc'][pos]
            else:
                merged_barcode += r2['pos2nuc'][pos]

        print(read.query_name, r1_barcode, r2_barcode, merged_barcode, sep='\t', file=tsv_out)

    bam.close()
    tsv_out.close()

def map_barcodes(intsv, oligo_info, outtsv):

    # Get barcode mismatches
    mismatch2barcode = dict()
    with open(oligo_info) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            line = line.strip().split('\t')
            line = dict(zip(header, line))
            for e1 in misc.edit1(line['barcode']):
                if e1 in mismatch2barcode and mismatch2barcode[e1] != line['barcode']:
                    raise Exception("Barcode mismatches overlap!")
                mismatch2barcode[e1] = line['barcode']

    # Map barcodes
    outtsv = open(outtsv, 'w')
    with open(intsv) as f:

        header = f.readline().strip().split('\t')
        print(header[0], 'mapped_barcode', sep='\t', file=outtsv)
        for line in f:
            line = line.strip().split('\t')
            line = dict(zip(header, line))

            if line['merged_barcode'] in mismatch2barcode:
                print(line['read_id'], mismatch2barcode[line['merged_barcode']], sep='\t', file=outtsv)
            elif line['r1_barcode'] in mismatch2barcode:
                print(line['read_id'], mismatch2barcode[line['r1_barcode']], sep='\t', file=outtsv)
            elif line['r2_barcode'] in mismatch2barcode:
                print(line['read_id'], mismatch2barcode[line['r2_barcode']], sep='\t', file=outtsv)

    outtsv.close()

if __name__ == '__main__':
    extract_barcodes(sys.argv[1], sys.argv[2], sys.argv[3])
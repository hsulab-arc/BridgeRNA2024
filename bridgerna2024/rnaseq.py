import subprocess
from bridgerna2024 import *
import pysam
from collections import defaultdict
from Bio import SeqIO

def _run_rnaseq_snakemake(workdir, threads):
    """Run the RNA-seq snakemake pipeline."""
    subprocess.run(['snakemake', '--snakefile', RNASEQ_SNAKEFILE_PATH, '--configfile',
                    RNASEQ_SNAKEMAKE_CONFIG_PATH, '-j', str(threads),
                    '--config',  'workdir={}'.format(workdir)], check=True)

def map_alignment_positions(aln, outtsv):

    sequences = dict()
    for rec in SeqIO.parse(aln, 'fasta'):
        sequences[rec.id] = str(rec.seq)

    with open(outtsv, 'w') as outfile:
        print('plasmid', 'orig_pos', 'aln_pos', 'nuc', sep='\t', file=outfile)
        for p in sequences:
            origpos = 0
            for alnpos,nuc in enumerate(sequences[p]):
                if nuc == '-':
                    continue
                print(p, origpos, alnpos, nuc, sep='\t', file=outfile)
                origpos += 1

def calculate_coverage(merged_bamfile, unmerged_bamfile, plasmid_fasta, annot_bed, depth_output, abundance_output):

    plasmid_seq = None
    for rec in SeqIO.parse(plasmid_fasta, 'fasta'):
        plasmid_seq = str(rec.seq)
        break

    annots = dict()
    with open(annot_bed) as infile:
        for line in infile:
            line = line.strip().split('\t')
            annots[line[3]] = (int(line[1]), int(line[2]))

    merged_bam = pysam.AlignmentFile(merged_bamfile, 'rb')
    ref_length = merged_bam.header['SQ'][0]['LN']
    fwd_coverage_unfiltered = dict()
    rev_coverage_unfiltered = dict()
    fwd_coverage_annot = dict()
    rev_coverage_annot = dict()
    for i in range(ref_length):
        fwd_coverage_unfiltered[i] = 0
        rev_coverage_unfiltered[i] = 0
        fwd_coverage_annot[i] = 0
        rev_coverage_annot[i] = 0

    for read in merged_bam:
        if read.is_unmapped:
            continue
        if read.is_reverse:
            for i in range(read.reference_start, read.reference_end):
                rev_coverage_unfiltered[i] += 1

            refstart_annot = False
            refend_annot = False
            for a in annots:
                if read.reference_start >= annots[a][0]:
                    refstart_annot = True
                if read.reference_end <= annots[a][1]:
                    refend_annot = True

            if refstart_annot and refend_annot:
                for i in range(read.reference_start, read.reference_end):
                    rev_coverage_annot[i] += 1
        else:
            for i in range(read.reference_start, read.reference_end):
                fwd_coverage_unfiltered[i] += 1

            refstart_annot = False
            refend_annot = False
            for a in annots:
                if read.reference_start >= annots[a][0]:
                    refstart_annot = True
                if read.reference_end <= annots[a][1]:
                    refend_annot = True

            if refstart_annot and refend_annot:
                for i in range(read.reference_start, read.reference_end):
                    fwd_coverage_annot[i] += 1

    merged_bam.close()

    unmerged_bam = pysam.AlignmentFile(unmerged_bamfile, 'rb')
    for read in unmerged_bam:
        if read.is_unmapped:
            continue
        if not read.is_read1:
            continue
        if read.mate_is_unmapped:
            continue
        tlen =  read.template_length
        if read.is_reverse:
            refstart = read.reference_end + tlen
            refend = read.reference_end
            for i in range(refstart, refend):
                rev_coverage_unfiltered[i] += 1

            refstart_annot = False
            refend_annot = False
            for a in annots:
                if refstart >= annots[a][0]:
                    refstart_annot = True
                if refend <= annots[a][1]:
                    refend_annot = True

            if refstart_annot and refend_annot:
                for i in range(read.reference_start, read.reference_end):
                    rev_coverage_annot[i] += 1

        else:
            refstart = read.reference_start
            refend = read.reference_start + tlen
            for i in range(refstart, refend):
                fwd_coverage_unfiltered[i] += 1

            refstart_annot = False
            refend_annot = False
            for a in annots:
                if refstart >= annots[a][0]:
                    refstart_annot = True
                if refend <= annots[a][1]:
                    refend_annot = True

            if refstart_annot and refend_annot:
                for i in range(read.reference_start, read.reference_end):
                    fwd_coverage_annot[i] += 1
    unmerged_bam.close()

    with open(depth_output, 'w') as f:
        print('pos', 'annot', 'fwd_unfilt', 'rev_unfilt', 'fwd_annot', 'rev_annot', sep='\t', file=f)
        for i in fwd_coverage_unfiltered:
            this_annot = '.'
            for a in annots:
                if i >= annots[a][0] and i <= annots[a][1]:
                    this_annot = a
            print(i, this_annot, fwd_coverage_unfiltered[i], rev_coverage_unfiltered[i], fwd_coverage_annot[i],
                  rev_coverage_annot[i], sep='\t', file=f)

    read_counts = defaultdict(int)
    merged_bam = pysam.AlignmentFile(merged_bamfile, 'rb')
    for read in merged_bam:
        if read.is_unmapped:
            continue
        strand = "+"
        if read.is_reverse:
            strand = '-'
        read_counts[(read.reference_start, read.reference_end, strand, read.get_reference_sequence())] += 1
    merged_bam.close()

    unmerged_bam = pysam.AlignmentFile(unmerged_bamfile, 'rb')
    for read in unmerged_bam:
        if read.is_unmapped:
            continue
        if not read.is_read1:
            continue
        if read.mate_is_unmapped:
            continue

        if read.is_reverse:
            strand = '-'
            refstart = read.reference_end + read.template_length
            refend = read.reference_end
        else:
            strand = '+'
            refstart = read.reference_start
            refend = read.reference_start + read.template_length

        if refstart == refend:
            print(read)

        read_counts[(refstart, refend, strand, plasmid_seq[refstart:refend])] += 1
    unmerged_bam.close()

    sorted_read_counts = list(reversed(sorted(read_counts.items(), key=lambda x: x[1])))

    aggregated_read_counts = []
    while True:
        if len(sorted_read_counts) == 0:
            break

        new_sorted_read_counts = []
        this_read = sorted_read_counts[0][0]
        this_count = sorted_read_counts[0][1]
        this_minstart = this_read[0]
        this_maxend = this_read[1]
        for rd, ct in sorted_read_counts[1:]:
            if rd[0] >= this_read[0] - 5 and rd[0] <= this_read[0] + 5 and rd[1] >= this_read[1] - 5 and rd[1] <= this_read[1] + 5 and rd[2] == this_read[2]:
                this_count += ct
                if rd[0] < this_minstart:
                    this_minstart = rd[0]
                if rd[1] > this_maxend:
                    this_maxend = rd[1]
            else:
                new_sorted_read_counts.append((rd, ct))

        out = {'start': this_read[0], 'end': this_read[1], 'strand': this_read[2], 'seq': this_read[3],
               'minstart': this_minstart, 'maxend': this_maxend, 'count': this_count}
        aggregated_read_counts.append(out)
        sorted_read_counts = list(new_sorted_read_counts)

    aggregated_read_counts = list(reversed(sorted(aggregated_read_counts, key=lambda x: x['count'])))

    with open(abundance_output, 'w') as f:
        print('start', 'end', 'strand', 'minstart', 'maxend', 'count', 'seq', sep='\t', file=f)
        for res in aggregated_read_counts:
            print(res['start'], res['end'], res['strand'], res['minstart'], res['maxend'], res['count'], res['seq'],
                  sep='\t', file=f)

if __name__ == '__main__':
    calculate_coverage(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

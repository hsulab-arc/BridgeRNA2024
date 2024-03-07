import gzip, re
from types import ModuleType
from hashlib import md5 as hashlib_md5
from base64 import b64encode

def hash_string(string):
    return b64encode(hashlib_md5(string.encode('utf-8')).digest()).decode('utf-8').rstrip('=')

def edit_distance(s1, s2):
    m, n = len(s1), len(s2)
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0:
                dp[i][j] = j
            elif j == 0:
                dp[i][j] = i
            elif s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]
            else:
                dp[i][j] = 1 + min(dp[i - 1][j], dp[i][j - 1], dp[i - 1][j - 1])

    return dp[m][n]

def edit1(s):
    nucs = ['A', 'C', 'G', 'T']
    for i in range(len(s)):
        for nuc in nucs:
            slist = list(s)
            slist[i] = nuc
            yield ''.join(slist)

def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_perc_identity(read):
    total = read.query_alignment_length
    identity = (total - read.get_tag('NM')) / total
    return identity

def fastq_reader(fastq_path):
    with gzip.open(fastq_path, 'rt') as f:
        while True:
            name = f.readline().rstrip()
            if name == '':
                break
            seq = f.readline().rstrip()
            spacer = f.readline().rstrip()
            qual = f.readline().rstrip()
            yield name, seq, spacer, qual

def batch_iterator(iterator, batch_size):
    batch = []
    for item in iterator:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

def count_reads(fastq):
    count = 0
    with gzip.open(fastq, 'rt') as f:
        for line in f:
            count += 1
    return count / 4

def count_lines(f, has_header=False):
    count = 0
    with open(f, 'rt') as f:
        if has_header:
            f.readline()
        for line in f:
            count += 1
    return count

def print_variables_dictionary(variables, ignore_dunder=True, ignore_module=True):
    for var in variables:
        if ignore_dunder and var.startswith("__"):
            continue
        if ignore_module and isinstance(variables[var], ModuleType):
            continue
        print(f"\t{var}: {variables[var]}")

def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

def retrieve_read_pairs_from_name_sorted_bam(sam):

    qry_name = None
    read1 = []
    read2 = []
    for read in sam:
        if qry_name is None:
            qry_name = read.query_name
        if read.query_name != qry_name:
            yield read1, read2
            qry_name = read.query_name
            read1 = []
            read2 = []

        if read.is_read1:
            read1.append(read)
        else:
            read2.append(read)
    yield read1, read2

def get_deletion_lengths(cigarstring):
    matches = re.findall(r'\d+D', cigarstring)
    for m in matches:
        yield int(m[:-1])

def get_insertion_lengths(cigarstring):
    matches = re.findall(r'\d+I', cigarstring)
    for m in matches:
        yield int(m[:-1])


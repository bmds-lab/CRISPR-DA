'''
utils.py

This file holder various utility methods that are used throughout the project.
'''
from datetime import datetime

def parse_fna(stream):
    '''Parse some iterable object as a multi-FASTA file.
    Yield after reading each FASTA block.

    Arguments:
        stream (iterable):  An iterable object to read
    '''
    header = None
    seqs = []
    for line in stream:
        line = line.strip()

        if line[0] == '>':
            if header is not None:
                yield header, ''.join(seqs)
            header = line
        else:
            seqs.append(line)
    yield header, ''.join(seqs)

def rc(dna):
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

# Function that replaces U with T in the sequence (to go back from RNA to DNA)
def trans_to_dna(rna: str):
    switch_UT = str.maketrans('U', 'T')
    dna = rna.translate(switch_UT)
    return dna

def one_hot_encode(seq, z='ATCG'):
    return [list(map(lambda x: 1.0 if x==c else 0.0, z)) for c in seq]

def colour_code(score: float):
    if score < 25:
        return '#e06666'
    elif score < 50:
        return '#f6b26b'
    elif score < 75:
        return '#ffe599'
    else:
        return '#93c47d'

def letter_code(score: float):
    if score < 25:
        return 'E'
    elif score < 50:
        return 'H'
    elif score < 75:
        return 'M'
    else:
        return 'L'

# Function that formats provided text with time stamp
def printer(stringFormat):
    print('>>> {}:\t{}\n'.format(
        datetime.now().strftime("%Y-%m-%d %H:%M:%S:%f"),
        stringFormat
    ))

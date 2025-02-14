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
    return [list(map(lambda x: 1 if x==c else 0, z)) for c in seq]

assert(one_hot_encode('ATCG') == [[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]])

def one_hot_decode(onehot):
    d = {
        0 : "A",
        1 : "T",
        2 : "C",
        3 : "G"
    }
    seq = ""
    for encoding in onehot:
        seq += d[encoding.index(1)]
    return seq

assert(one_hot_decode([[1,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,1]]) == 'ATCG')
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
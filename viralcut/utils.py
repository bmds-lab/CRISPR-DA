

def rc(dna):
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

# Function that replaces U with T in the sequence (to go back from RNA to DNA)
def trans_to_dna(rna: str):
    switch_UT = str.maketrans('U', 'T')
    dna = rna.translate(switch_UT)
    return dna
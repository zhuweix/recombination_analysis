def rc_seq(seq=""):
    """Returns the reverse compliment sequence."""
    rc_dict = {
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c",
        "n": "n",
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N"
    }
    rc_nt_ls = map(lambda x: rc_dict[x], reversed(seq))
    rc_seq_ = "".join(rc_nt_ls)
    return rc_seq_

def translate_exon(seq):
    """Translate exon by normal codon."""

    codon = {
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "CTT": "L",
        "CTC": "L",
        "CTG": "L",
        "CTA": "L",
        "TTA": "L",
        "TTG": "L",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "TTT": "F",
        "TTC": "F",
        "ATG": "M",
        "TGT": "C",
        "TGC": "C",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACT": "T",
        "ACC": "T",
        "ACG": "T",
        "ACA": "T",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "AGT": "S",
        "AGC": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TGG": "W",
        "CAA": "Q",
        "CAG": "Q",
        "AAT": "N",
        "AAC": "N",
        "CAT": "H",
        "CAC": "H",
        "GAA": "E",
        "GAG": "E",
        "GAT": "D",
        "GAC": "D",
        "AAA": "K",
        "AAG": "K",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGA": "R",
        "AGG": "R",
        "TAA": "*",
        "TAG": "*",
        "TGA": "*"
    }

    seq = seq.upper()
    size = len(seq)
    if size % 3 != 0:
        print("Length is {}, invalid !".format(size))
        return ''
    protein = [codon[seq[id_: id_+3]] for id_ in range(0, size, 3)]
    protein = "".join(protein)
    return protein
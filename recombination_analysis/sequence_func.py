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


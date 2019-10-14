import os

def simple_fasta_write(fname, names, seqs, linewidth=80):
    """Write the fasta file."""
    with open(fname, "w") as fasta_p:
        total = len(names)
        for _id in range(total):
            fasta_p.write(">%s\n" % names[_id])
            _seq = seqs[_id]
            # seq with linebreak
            _seq_wb = "\n".join([_seq[i:linewidth+i]
                                 for i in range(0, len(_seq), linewidth)])
            fasta_p.write(_seq_wb)
            fasta_p.write("\n")


def simple_fastq_load(fname: str) -> tuple:
    """Load Fastq file with Phred Score in 32-ASCII code."""
    names = []
    seqs = []
    scores = []
    if not os.path.isfile(fname):
        print("%s is not available!" % fname)
        return None
    with open(fname) as filep:
        _cur_name = ""
        _cur_seq = []
        _cur_score = []
        for count, line in enumerate(filep):
            line = line.strip()
            if not line:
                continue
            # New entry
            if count % 4 == 0:
                if line[0] != '@':
                    print("Fastq file reading error in reading %s" %
                             line[:20])
                    return None
                if _cur_name:
                    names.append(_cur_name)
                    seqs.append(_cur_seq)
                    scores.append(_cur_score)
                _cur_name = line[1:]
                _cur_seq = []
                _cur_score = []
            elif count % 4 == 1:
                _cur_seq = line
            # load Phred Score
            elif count % 4 == 3:
                _cur_score = [ord(phred) - 33 for phred in line]
                if len(_cur_score) != len(_cur_seq):
                    print(
                        "Length of sequence and score doesnt match for %s" % _cur_name)
                    return None
        names.append(_cur_name)
        seqs.append(_cur_seq)
        scores.append(_cur_score)
    fq_dict = {}
    for n, seq, sc in zip(names, seqs, scores):
        fq_dict[n] = {'seq': seq, 'score': sc}
    return fq_dict


def simple_fasta_load(fname: str) -> tuple:
    """Load fasta file.
    Args:
    fname: name of the fasta file
    Returns:
    dict = {name: seq}
    """

    names = []
    seqs = []
    if not os.path.isfile(fname):
        print("%s is not available!" % fname)
        return None

    with open(fname, 'r') as fasta_p:
        _cur_name = ""
        _cur_seq = []
        for line in fasta_p:
            line_ = line.strip()
            # Skip empty lines
            if not line_:
                continue
            # New entry
            if line_[0] == ">":
                # Add previous entry
                if _cur_name:
                    names.append(_cur_name)
                    seqs.append("".join(_cur_seq))
                _cur_name = line_[1:]  # Omit >
                _cur_seq = []
            else:
                # Update sequence of the entry
                if not _cur_name:
                    print("One seq without entry")
                    print(line_)
                    return None
                _cur_seq.append(line_)

        # Update the final entry
        if _cur_name:
            names.append(_cur_name)
            seqs.append("".join(_cur_seq))
    fa_dict = {n: s for n, s in zip(names, seqs)}
    return fa_dict


def simple_gff3_load(fname, return_fasta=False):
    """Load gff3 files."""
    entries = {}
    with open(fname, "r") as gff3:
        for line in gff3:
            line = line.strip()
            if line[0] == "#":
                if line == "##FASTA":
                    break
                continue
            entry = line.split("\t")
            if len(entry) < 9:
                print(
                    "Error loading %s: Less than 9 items in the entry\n%s" % (fname, line))
                continue
            chrom = entry[0]
            author = entry[1]
            type_ = entry[2]
            start = int(entry[3]) - 1
            end = int(entry[4])
            direction = entry[6]
            note = {}
            note_entries = entry.pop().split(";")
            for n in note_entries:
                class_ = n.split('=')[0]
                content = n.split('=')[1]
                note[class_] = content
            entries.setdefault(chrom, [])
            entries[chrom].append({
                'type': type_,
                "start": start,
                'end': end,
                'direction': direction,
                'note': note,
                'author': author})
        print("%s: Gff entries are analyzed" % fname)

    if not return_fasta:
        return entries
    else:
        names = []
        seqs = []
        with open(fname, "r") as gff3:
            line = gff3.readline()
            line = line.strip()
            while not line == "##FASTA":
                line = gff3.readline()
                line = line.strip()
            _cur_name = ""
            _cur_seq = []
            line = gff3.readline()
            while line:
                line_ = line.strip()
                line = gff3.readline()
                # Skip empty lines
                if not line_:
                    continue
                # New entry
                if line_[0] == ">":
                    # Add previous entry
                    if _cur_name:
                        names.append(_cur_name)
                        seqs.append("".join(_cur_seq))
                    _cur_name = line_[1:]  # Omit >
                    _cur_seq = []
                else:
                    # Update sequence of the entry
                    if not _cur_name:
                        print("One seq without entry")
                        print(line_)
                        return entries, None
                    _cur_seq.append(line_)
            # Update the final entry
            if _cur_name:
                names.append(_cur_name)
                seqs.append("".join(_cur_seq))
        fa_dict = {n: s for n, s in zip(names, seqs)}
        return entries, fa_dict


def simple_bed_load(fname):
    """Load bed files."""
    entries = {}
    with open(fname, "r") as gff3:
        for line in gff3:
            line = line.strip()
            if line.startswith('track'):
                continue
            elif line.startswith('#'):
                continue
            elif not line:
                continue
            entry = line.split("\t")
            if len(entry) < 3:
                print(
                    "Error loading %s: Less than 3 items in the entry\n%s" % (fname, line))
                continue
            chrom = entry[0]
            start = int(entry[1])
            end = int(entry[2])
            possible_extra_columns = ['name', 'score', 'direction']
            cur_entry = {'start': start, 'end': end}
            for i in range(min(len(possible_extra_columns), len(entry) - 3)):
                cur_entry[possible_extra_columns[i]] = entry[3 + i]

            entries.setdefault(chrom, [])
            entries[chrom].append(cur_entry)

        print("%s: bed entries are analyzed" % fname)
    return entries


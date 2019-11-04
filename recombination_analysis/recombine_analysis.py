from .biofile_func import simple_fasta_load
import os
import sys
from bitarray import bitarray as bbitarray
from matplotlib import cm
from matplotlib import colors as mcolors
from matplotlib import colorbar as mcolorbar
from matplotlib import pyplot as plt
import numpy as np


def gen_ref_gene_dict(fn: str, window=20):
    genome = simple_fasta_load(fn)
    gene_id = {}
    total = len(genome)
    rkmer_dict = {}

    for i, g in enumerate(sorted(genome.items())):
        name, seq = g
        seq = seq.upper()
        gene_id[i] = name

        for p in range(len(seq) - window + 1):
            kmer = seq[p: p + window]
            if kmer not in rkmer_dict:
                gene_no = total * bbitarray('0')
                gene_no[i] = True
                rkmer_dict[kmer] = gene_no
            else:
                rkmer_dict[kmer][i] = True
    return rkmer_dict, gene_id


def gen_ref_gene_dict_inter(genome: dict, window=20):
    gene_id = {}
    total = len(genome)
    rkmer_dict = {}

    for i, g in enumerate(sorted(genome.items())):
        name, seq = g
        seq = seq.upper()
        gene_id[i] = name

        for p in range(len(seq) - window + 1):
            kmer = seq[p: p + window]
            if kmer not in rkmer_dict:
                gene_no = total * bbitarray('0')
                gene_no[i] = True
                rkmer_dict[kmer] = gene_no
            else:
                rkmer_dict[kmer][i] = True
    return rkmer_dict, gene_id


def bidirect_assign_path(seq: str, rkmer_dict: dict, window: int, step: int):
    assign_path_forward = []
    last_i = 0
    assign_path_backward = []
    first_i = 0
    for i in range(0, len(seq) - window + 1, step):
        kmer = seq[i: i + window]
        if kmer not in rkmer_dict:
            continue
        last_i = i
        assign_no = rkmer_dict[kmer]
        if not assign_path_forward:
            assign_path_forward = [
                len(assign_no) * bbitarray('0'), [i, assign_no]]
            continue
        orth = assign_path_forward[-1][1] & assign_no
        if any(orth):
            assign_path_forward[-1][1] = orth
        else:
            assign_path_forward.append([i, assign_no])

    for p in assign_path_forward[1:]:
        assign_path_forward[0] |= p[1]
    assign_path_forward.append((last_i, None))

    for i in range(len(seq) - window + 1, -1,  -step):
        kmer = seq[i: i + window]
        if kmer not in rkmer_dict:
            continue
        first_i = i
        assign_no = rkmer_dict[kmer]
        if not assign_path_backward:
            assign_path_backward = [
                len(assign_no) * bbitarray('0'), [i, assign_no]]
            continue
        orth = assign_path_backward[-1][1] & assign_no
        if any(orth):
            assign_path_backward[-1][1] = orth
        else:
            assign_path_backward.append([i, assign_no])

    for p in assign_path_backward[1:]:
        assign_path_backward[0] |= p[1]
    assign_path_backward.append((first_i, None))
    # adjust backward sequence
    assign_path_backward = [assign_path_backward[0]] + \
        assign_path_backward[-1:0:-1]
    return assign_path_forward, assign_path_backward


def merge_direct_assign(path_for, path_back):
    # Intersect == empty  => homology region
    # Intersect != empty  => recombined region
    assigned_genes = len(path_for[0]) * bbitarray('0')
    aleft = max(path_for[1][0], path_back[1][0])
    aright = min(path_for[-1][0], path_back[-1][0])
    assign_path = [assigned_genes]
    cur_left = -1
    cur_right = -1
    cur_assign = None
    cur_fid = 1
    cur_bid = 1
    for p in range(aleft, aright):
        is_update = False
        # Locate forward assign
        cur_f1 = path_for[cur_fid]
        cur_f2 = path_for[cur_fid + 1]
        if p >= cur_f2[0]:
            cur_fid += 1
            cur_f1 = path_for[cur_fid]
            cur_f2 = path_for[cur_fid + 1]
            is_update = True
        cur_for_assign = cur_f1[1]
        # Locate backward assign
        cur_b1 = path_back[cur_bid]
        cur_b2 = path_back[cur_bid + 1]
        if p >= cur_b2[0]:
            cur_bid += 1
            cur_b1 = path_back[cur_bid]
            cur_b2 = path_back[cur_bid + 1]
            is_update = True
        cur_back_assign = cur_b2[1]
        if (cur_assign is None) or is_update:
            intersect = cur_for_assign & cur_back_assign
            if cur_left >= 0:
                assign_path.append([cur_left, cur_right, cur_assign])
            cur_left = p
            cur_right = p
            if intersect.any():
                cur_assign = [intersect, None]
                assign_path[0] |= intersect
            else:
                cur_assign = [cur_for_assign, cur_back_assign]
                assign_path[0] |= (cur_for_assign | cur_back_assign)
        else:
            cur_right = p
    assign_path.append([cur_left, cur_right, cur_assign])

    return assign_path


def count_bidirect_kmer_evidence(seq, rkmer_dict, window, assign_path):
    # Count forward, backward, and overlapped evidence
    share_cur_count = []
    share_cur_ratio = []

    for ap in assign_path[1:]:
        cur_count = 0

        for p in range(ap[0], ap[1] + 1):
            kmer = seq[p: p + window]
            # current
            if kmer in rkmer_dict:
                cur_count += 1
        share_cur_count.append(cur_count)
        share_cur_ratio.append(cur_count / (ap[1] - ap[0] + 1))

    return share_cur_ratio, share_cur_count


def determine_recombine_bidirect(seq: str, rkmer_dict: dict, window: int, step: int):
    seq = seq.upper()
    ap_forward, ap_backward = bidirect_assign_path(
        seq, rkmer_dict, window, step)
    if not ap_forward:
        return None, None, None
    if ap_forward == [(0, None)]:
        return None, None, None
    ap = merge_direct_assign(ap_forward, ap_backward)

    share_cur_ratio, share_cur_count = count_bidirect_kmer_evidence(
        seq, rkmer_dict, window, ap)

    return ap, share_cur_ratio, share_cur_count


def draw_bidirect_gene_assign(assign_path, rgene_id, gene, length,
                              figname, share_cur_ratio):
    cmap = cm.get_cmap('tab10')
    fig = plt.figure(figsize=(6.5, 6.5), dpi=300)

    text_pad = 10
    ax = plt.subplot2grid((12, 1), (0, 0))
    cbar = mcolorbar.ColorbarBase(ax, cmap=cmap, orientation='horizontal',
                                  ticks=np.linspace(0, 1, 11),
                                  ticklocation='top')
    ax = plt.subplot2grid((12, 1), (1, 0), rowspan=11)

    plt.sca(ax)
    plt.xlim(0, length)

    assigned_ref_ids = np.nonzero(assign_path[0])[0]
    assigned_ref_ids = {id_: i for i, id_ in enumerate(assigned_ref_ids)}
    assigned_ref_gene = {i: rgene_id[id_]
                         for id_, i in assigned_ref_ids.items()}
    total = len(assigned_ref_ids)

    for i, p in enumerate(assign_path[1:]):
        cur_start, cur_end, cur_assign = p
        c = share_cur_ratio[i]
        cur_id1 = [assigned_ref_ids[id_]
                   for id_ in np.nonzero(cur_assign[0])[0]]
        cur_id2 = None
        if cur_assign[1]:
            cur_id2 = [assigned_ref_ids[id_]
                       for id_ in np.nonzero(cur_assign[1])[0]]

        if not cur_id2:
            for id_ in cur_id1:
                ax.plot((cur_start, cur_end), (id_ + .5,
                                               id_ + .5), lw=1.5, color=cmap(c))
        else:
            for id_ in cur_id1:
                ax.plot((cur_start, cur_end), (id_ + .5, id_ + .5),
                        lw=1.5, ls=':', color=cmap(c))
            for id_ in cur_id2:
                ax.plot((cur_start, cur_end), (id_ + .5, id_ + .5),
                        lw=1.5, ls=':', color=cmap(c))

        if i != len(assign_path) - 2:
            c = share_cur_ratio[i]
            ax.plot((cur_end, cur_end), (0, len(assigned_ref_ids)), color=cmap(c),
                    lw=1)
    for i in range(total):
        plt.text(length + text_pad, i + .5,
                 assigned_ref_gene[i], ha='left', va='center')
    plt.ylim([0, total])
    plt.xlabel('%s ORF position (nt)' % gene)
    plt.gca().get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.savefig('%s.png' % figname)
    plt.close()


def draw_gene_assign_dotplot(assign_path, rgene_id: dict, gene: str,
                             figname: str, qseq: str, rseqs: dict, oseq: str, oname: str,
                             min_rec: int, top=0, ksize=30):

    assigned_ref_ids = np.nonzero(assign_path[0])[0]
    assigned_ref_ids = {id_: i for i, id_ in enumerate(assigned_ref_ids)}
    assigned_ref_gene = {i: rgene_id[id_]
                         for id_, i in assigned_ref_ids.items()}
    if gene not in assigned_ref_gene.values():
        return
    ref_locations = {}
    ref_rec_length = {}
    ref_homo_locations = {}
    self_locations = []
    rec_locations = []
    homo_locations = []
    self_id = [i for i, name in assigned_ref_gene.items() if name == gene][0]
    total = len(assigned_ref_ids)
    length = len(qseq)

    for i, p in enumerate(assign_path[1:]):
        cur_start, cur_end, cur_assign = p
        cur_id1 = [assigned_ref_ids[id_]
                   for id_ in np.nonzero(cur_assign[0])[0]]
        cur_id2 = None
        if cur_assign[1]:
            cur_id2 = [assigned_ref_ids[id_]
                       for id_ in np.nonzero(cur_assign[1])[0]]

        if not cur_id2:
            if self_id in cur_id1:
                self_locations.append((cur_start, cur_end))
            else:
                if cur_end - cur_start >= min_rec:
                    rec_locations.append((cur_start, cur_end))
                    for id_ in cur_id1:
                        ref_locations.setdefault(id_, [])
                        ref_locations[id_].append((cur_start, cur_end))
                        ref_rec_length.setdefault(id_, 0)
                        ref_rec_length[id_] += (cur_end - cur_start)
        else:
            homo_locations.append((cur_start, cur_end))
            all_ids = [id_ for id_ in cur_id1]
            all_ids.extend(cur_id2)
            all_ids = list(set(all_ids))
            for id_ in all_ids:
                if id_ != self_id:
                    ref_homo_locations.setdefault(id_, [])
                    ref_homo_locations[id_].append((cur_start, cur_end))

    ref_rec_length = [(id_, length) for id_, length in ref_rec_length.items()]
    ref_rec_length.sort(key=lambda x: x[1], reverse=True)
    if top != 0:
        ref_rec_length = ref_rec_length[:top]
    if not ref_rec_length:
        return

    # Draw self-matrix
    fig = plt.figure(figsize=(6.5, 6.5), dpi=300)
    plt.xlim(0, length)
    plt.ylim(0, length)
    p = 0
    kdict = {}
    for i in range(len(qseq) - ksize + 1):
        kmer = qseq[i: i+ksize]
        if kmer not in kdict:
            kdict[kmer] = [i]
        else:
            kdict[kmer].append(i)
    self_x = []
    self_y = []
    rec_x = []
    rec_y = []
    homo_x = []
    homo_y = []
    all_x = []
    all_y = []
    for locs in kdict.values():
        sloc = []
        rloc = []
        hloc = []
        for reg in self_locations:
            sloc.extend([l for l in locs if reg[0] <= l <= reg[1]])
        for reg in rec_locations:
            rloc.extend([l for l in locs if reg[0] <= l <= reg[1]])
        for reg in homo_locations:
            hloc.extend([l for l in locs if reg[0] < l < reg[1]])
        self_x.extend(sloc * len(sloc))
        self_y.extend([s for s in sloc for t in sloc])
        rec_x.extend(rloc * len(rloc))
        rec_y.extend([s for s in rloc for t in rloc])
        homo_x.extend(hloc * len(hloc))
        homo_y.extend([s for s in hloc for t in hloc])
        all_x.extend(locs * len(locs))
        all_y.extend([u for u in locs for v in locs])
    plt.scatter(all_x, all_y, s=.5, color='lightsteelblue')
    plt.scatter(homo_x, homo_y, color='cyan', s=.5)
    plt.scatter(self_x, self_y, color='blue', s=.5, label='Non-recombined')
    plt.scatter(rec_x, rec_y, color='orange', s=.5)
    plt.xlabel('{} (nt)'.format(gene))
    plt.ylabel('{} (nt)'.format(gene))
    lrec = mpatches.Patch(color='orange', label='Recombined region')
    lhomo = mpatches.Patch(color='cyan', label='Homologous flanking region')
    lnon = mpatches.Patch(color='blue', label='Non-recombined region')
    lo = mpatches.Patch(color='lightsteelblue', label='Other')
    plt.legend(handles=[lrec, lhomo, lnon, lo], loc=2, fontsize=9)
    plt.tight_layout()
    plt.savefig('%s.0.png' % figname)
    plt.close()

    # Draw others
    # Generate Kmer from ortholog
    okmer_dict = {}
    for i in range(len(oseq) - ksize + 1):
        kmer = oseq[i: i + ksize]
        okmer_dict.setdefault(kmer, [])
        okmer_dict[kmer].append(i)

    for pid, (rid, _) in enumerate(ref_rec_length, 1):
        # recombined kmers
        rec_kmer_dict = {}
        rec_anchor_dict = {}
        rname = assigned_ref_gene[rid]
        rseq = rseqs[rname]
        ref_rec = ref_locations[rid]

        ref_dict = {}
        # ref kmer dict
        for i in range(len(rseq) - ksize + 1):
            kmer = rseq[i: i + ksize]
            ref_dict.setdefault(kmer, [])
            ref_dict[kmer].append(i)

        # Obtain kmers from query
        # recombined kmer
        for reg in ref_rec:
            for i in range(reg[0], reg[1] + 1):
                kmer = qseq[i: i + ksize]
                if kmer in ref_dict:
                    if kmer not in okmer_dict:
                        rec_anchor_dict.setdefault(kmer, [])
                        rec_anchor_dict[kmer].append(i)
                    else:
                        rec_kmer_dict.setdefault(kmer, [])
                        rec_kmer_dict[kmer].append(i)

        fig = plt.figure(figsize=(6.5, 6.5), dpi=300)
        # Ref-query plot
        max_size = max(len(rseq), len(qseq), len(oseq))
        plt.xlim(0, max_size)
        plt.ylim(0, max_size)
        x_pos = []
        y_pos = []
        rec_x = []
        rec_y = []
        an_x = []
        an_y = []
        for kmer in kdict:
            if kmer in rec_anchor_dict:
                qloc = rec_anchor_dict[kmer]
                rloc = ref_dict[kmer]
                an_x.extend(rloc * len(qloc))
                an_y.extend([p for p in qloc for q in rloc])
        for kmer in kdict:
            if kmer in rec_kmer_dict:
                qloc = rec_kmer_dict[kmer]
                rloc = ref_dict[kmer]
                rec_x.extend(rloc * len(qloc))
                rec_y.extend([p for p in qloc for q in rloc])
        for kmer in kdict:
            if kmer in ref_dict:
                qloc = kdict[kmer]
                rloc = ref_dict[kmer]
                x_pos.extend(rloc * len(qloc))
                y_pos.extend([p for p in qloc for q in rloc])
        plt.scatter(x_pos, y_pos, color='grey', s=.5)
        plt.scatter(rec_x, rec_y, color='orange', s=.5)
        plt.scatter(an_x, an_y, color='red', s=.5)
        plt.xlabel('{} (nt)'.format(rname))
        plt.ylabel('{} (nt)'.format(gene))
        lanchor = mpatches.Patch(color='red', label='Recombined kmer')
        lrec = mpatches.Patch(color='orange', label='Recombined region')
        lhomo = mpatches.Patch(color='grey', label='Shared region')
        plt.legend(handles=[lanchor, lrec, lhomo], loc=2, fontsize=9)
        plt.tight_layout()
        plt.savefig('{}.{}.{}.a.png'.format(figname, pid, rname))
        plt.close()

        # Query-ortholog
        fig = plt.figure(figsize=(6.5, 6.5), dpi=300)
        plt.xlim(0, max_size)
        plt.ylim(0, max_size)
        x_pos = []
        y_pos = []

        for kmer in kdict:
            if kmer in okmer_dict:
                qloc = kdict[kmer]
                rloc = okmer_dict[kmer]
                x_pos.extend(rloc * len(qloc))
                y_pos.extend([p for p in qloc for q in rloc])
        plt.scatter(x_pos, y_pos, color='grey', s=.5)
        lhomo = mpatches.Patch(color='grey', label='Shared region')
        plt.legend(handles=[lhomo], loc=2, fontsize=9)
        plt.xlabel('{} (nt)'.format(oname))
        plt.ylabel('{} (nt)'.format(gene))
        plt.tight_layout()
        plt.savefig('{}.{}.{}.b.png'.format(figname, pid, rname))
        plt.close()

        # Ref-ref
        fig = plt.figure(figsize=(6.5, 6.5), dpi=300)
        plt.xlim(0, max_size)
        plt.ylim(0, max_size)
        x_pos = []
        y_pos = []
        rec_x = []
        rec_y = []
        an_x = []
        an_y = []
        for kmer in ref_dict:
            if kmer in rec_anchor_dict:
                qloc = ref_dict[kmer]
                an_x.extend(qloc * len(qloc))
                an_y.extend([p for p in qloc for q in qloc])
            if kmer in rec_kmer_dict:
                qloc = ref_dict[kmer]
                rec_x.extend(qloc * len(qloc))
                rec_y.extend([p for p in qloc for q in qloc])
            else:
                qloc = ref_dict[kmer]
                x_pos.extend(qloc * len(qloc))
                y_pos.extend([p for p in qloc for q in qloc])

        plt.scatter(x_pos, y_pos, color='lightsteelblue', s=.5)
        plt.scatter(rec_x, rec_y, color='orange', s=.5)
        plt.scatter(an_x, an_y, color='red', s=1)
        plt.xlabel('{} (nt)'.format(rname))
        plt.ylabel('{} (nt)'.format(rname))

        lanchor = mpatches.Patch(color='red', label='Recombined kmer')
        lrec = mpatches.Patch(color='orange', label='Recombined region')
        lhomo = mpatches.Patch(color='lightsteelblue', label='Other')
        plt.legend(handles=[lanchor, lrec, lhomo], loc=2, fontsize=9)
        plt.tight_layout()
        plt.savefig('{}.{}.{}.c.png'.format(figname, pid, rname))
        plt.close()

from .biofile_func import simple_fasta_load
import os
import sys
from bitarray import bitarray
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
                gene_no = total * bitarray('0')
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
                gene_no = total * bitarray('0')
                gene_no[i] = True
                rkmer_dict[kmer] = gene_no
            else:
                rkmer_dict[kmer][i] = True
    return rkmer_dict, gene_id

def gen_assign_path(seq: str, rkmer_dict: dict, window: int, step: int):
    assign_path = []
    last_i = 0
    for i in range(0, len(seq) - window + 1, step):
        kmer = seq[i: i + window]
        if kmer not in rkmer_dict:
            continue
        last_i = i
        assign_no = rkmer_dict[kmer]
        if not assign_path:
            assign_path = [len(assign_no) * bitarray('0'), [i, assign_no]]
            continue
        orth = assign_path[-1][1] & assign_no
        if any(orth):
            assign_path[-1][1] = orth
        else:
            assign_path.append([i, assign_no])

    for p in assign_path[1:]:
        assign_path[0] |= p[1]
    assign_path.append((last_i, None))

    return assign_path

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
            assign_path_forward = [len(assign_no) * bitarray('0'), [i, assign_no]]
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
            assign_path_backward = [len(assign_no) * bitarray('0'), [i, assign_no]]
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
    assign_path_backward = [assign_path_backward[0]] + assign_path_backward[-1:0:-1]
    return assign_path_forward, assign_path_backward

def merge_direct_assign(path_for, path_back):
    # Intersect == empty  => homology region
    # Intersect != empty  => recombined region
    assigned_genes= len(path_for[0]) * bitarray('0')
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
            if cur_left >=0:
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
    share_cur_score = []
    share_pre_score = [None]
    seq_len = len(seq)
    # Count current
    ap = assign_path[1]
    gap_penalty = 0
    score = ap[1] - ap[0]
    for p in range(ap[0], ap[1] - window + 1):
        kmer = seq[p: p+window]
        if kmer in rkmer_dict:
            if gap_penalty >= window:
                score -= (gap_penalty - window + 1)
            gap_penalty = 0
        else:
            gap_penalty += 1
    if gap_penalty >= window:
        score -= (gap_penalty - window + 1)
    share_cur_score.append(score)
    if len(assign_path) == 2:
        return share_cur_score, share_pre_score
    else:
        pre_assign = ap[2]
        if pre_assign[1]:
            pre_assign = pre_assign[0] | pre_assign[1]
        else:
            pre_assign = pre_assign[0]
        for pid in range(2, len(assign_path)):
            ap = assign_path[pid]
            gap_penalty = 0
            score = ap[1] - ap[0]
            for p in range(ap[0], ap[1] - window + 1):
                kmer = seq[p: p+window]
                if kmer in rkmer_dict:
                    if gap_penalty >= window:
                        score -= (gap_penalty - window + 1)
                    gap_penalty = 0
                else:
                    gap_penalty += 1
            if gap_penalty >= window:
                score -= (gap_penalty - window + 1)                
            share_cur_score.append(score)
            gap_penalty = 0
            score = ap[1] - ap[0]
            for p in range(ap[0], ap[1] - window + 1):
                kmer = seq[p: p+window]
                if kmer in rkmer_dict:
                    if (rkmer_dict[kmer] & pre_assign).any():
                        if gap_penalty >= window:
                            score -= (gap_penalty - window + 1)
                        gap_penalty = 0
                    else:
                        gap_penalty += 1
                else:
                    gap_penalty += 1
            if gap_penalty >= window:
                score -= (gap_penalty - window + 1)
                gap_penalty = 0
            share_pre_score.append(score)
            pre_assign = ap[2]
            if pre_assign[1]:
                pre_assign = pre_assign[0] | pre_assign[1]
            else:
                pre_assign = pre_assign[0]            
    return share_cur_score, share_pre_score

def count_kmer_evidence(seq: str, rkmer_dict: dict, window: int, step: int, assign_path: list):
    # First count kmer in assigned pathes
    share_cur_kmer_count = []
    share_pre_kmer_count = [None]
    seq_len = len(seq)
    if len(assign_path) == 3:
        # count cur share
        share_cur_kmer_count.append(0)
        for i in range(0, len(seq) - window, step):
            kmer = seq[i:i+ window]
            if kmer not in rkmer_dict:
                continue
            share_cur_kmer_count[0] += 1
        return share_cur_kmer_count, share_pre_kmer_count
    for i in range(1, len(assign_path) -1):
        cur_count = 0
        pre_count = 0
        for j in range(assign_path[i][0], assign_path[i+1][0]):
            kmer = seq[j: j+window]
            if kmer not in rkmer_dict:
                continue
            cur_count += 1
            if i > 1:
                pre_assign = assign_path[i-1][1]
                pre_count += (pre_assign & rkmer_dict[kmer]).any()
        share_cur_kmer_count.append(cur_count)
        if i != 1:
            share_pre_kmer_count.append(pre_count)
        
    return share_cur_kmer_count, share_pre_kmer_count    
    

def determine_recombine(seq: str, rkmer_dict: dict, window: int, step: int):
    seq = seq.upper()
    assign_path = gen_assign_path(seq, rkmer_dict, window, step)
    if not assign_path:
        return None, None, None
    share_cur_kmer_count, share_pre_kmer_count = count_kmer_evidence(
        seq, rkmer_dict, window, step, assign_path)
    return assign_path, share_cur_kmer_count, share_pre_kmer_count

def determine_recombine_bidirect(seq: str, rkmer_dict: dict, window: int, step: int):
    seq = seq.upper()
    ap_forward, ap_backward = bidirect_assign_path(seq, rkmer_dict, window, step)
    if not ap_forward:
        return None, None, None
    if ap_forward == [(0, None)]:
        return None, None, None
    ap = merge_direct_assign(ap_forward, ap_backward)
    share_cur_score, share_pre_score = count_bidirect_kmer_evidence(
        seq, rkmer_dict, window, ap)
    return ap, share_cur_score, share_pre_score
    
    

def draw_gene_assign(assign_path, rgene_id: dict, gene: str, length: int,
                     figname: str, share_cur_kmer_count, share_pre_kmer_count):
    cmap = cm.get_cmap('tab10')
    fig = plt.figure(figsize=(8, 8), dpi=300)
    
    text_pad = 10
    ax = plt.subplot2grid((12, 1), (0,0))
    plt.title('%s' % gene)    
    cbar = mcolorbar.ColorbarBase(ax, cmap=cmap, orientation='horizontal', 
                                  ticks=np.linspace(0,1,11), 
                                  ticklocation='top') 
    ax = plt.subplot2grid((12, 1), (1,0), rowspan=11)

    plt.sca(ax)
    plt.xlim(0, length)

    assigned_ref_ids = np.nonzero(assign_path[0])[0]
    assigned_ref_ids = {id_: i for i, id_ in enumerate(assigned_ref_ids)}
    assigned_ref_gene = {i: rgene_id[id_]
                         for id_, i in assigned_ref_ids.items()}
    total = len(assigned_ref_ids)
    cur_start = assign_path[1][0]
    cur_id = [assigned_ref_ids[id_]
              for id_ in np.nonzero(assign_path[1][1])[0]]

    i = 0
    for i, p in enumerate(assign_path[2:]):
        cur_end = p[0]
        c =  min(share_cur_kmer_count[i] / (cur_end - cur_start+1), 1)        
        for id_ in cur_id:
            ax.plot((cur_start, cur_end), (id_ + .5, id_ + .5), lw=2, color=cmap(c))
        cur_start = cur_end
        cur_id = [assigned_ref_ids[id_]
                  for id_ in np.nonzero(p[1])[0]]
        if i != len(assign_path) -3:
            c = 1 - share_pre_kmer_count[i+1] / (share_cur_kmer_count[i+1] + 0.001)
            ax.plot((cur_start, cur_start), (0, len(assigned_ref_ids)), color=cmap(c),
                    lw=1)
        else:
            i += 1
    for i in range(total):
        plt.text(length + text_pad, i + .5,
                 assigned_ref_gene[i], ha='left', va='center')
    plt.xlabel('Gene orf position (nt)')
 
    plt.tight_layout()
    plt.savefig('%s.png' % figname)
    plt.close()

def draw_bidirect_gene_assign(assign_path, rgene_id, gene, length,
                              figname, share_cur_score, share_pre_score):
    cmap = cm.get_cmap('tab10')
    fig = plt.figure(figsize=(12, 12), dpi=72)
    
    text_pad = 10
    ax = plt.subplot2grid((12, 1), (0,0))
    plt.title('%s' % gene)    
    cbar = mcolorbar.ColorbarBase(ax, cmap=cmap, orientation='horizontal', 
                                  ticks=np.linspace(0,1,11), 
                                  ticklocation='top') 
    ax = plt.subplot2grid((12, 1), (1,0), rowspan=11)

    plt.sca(ax)
    plt.xlim(0, length)

    assigned_ref_ids = np.nonzero(assign_path[0])[0]
    assigned_ref_ids = {id_: i for i, id_ in enumerate(assigned_ref_ids)}
    assigned_ref_gene = {i: rgene_id[id_]
                         for id_, i in assigned_ref_ids.items()}
    total = len(assigned_ref_ids)

    for i, p in enumerate(assign_path[1:]):
        cur_start, cur_end, cur_assign = p
        c =  min(share_cur_score[i] / (cur_end - cur_start+1), 1)
        cur_id1 = [assigned_ref_ids[id_]
                  for id_ in np.nonzero(cur_assign[0])[0]]
        cur_id2 = None
        if cur_assign[1]:
            cur_id2 = [assigned_ref_ids[id_]
                      for id_ in np.nonzero(cur_assign[1])[0]]

        if not cur_id2:
            for id_ in cur_id1:
                ax.plot((cur_start, cur_end), (id_ + .5, id_ + .5), lw=1.5, color=cmap(c))
        else:
            for id_ in cur_id1:
                ax.plot((cur_start, cur_end), (id_ + .5, id_ + .5), lw=1.5, ls=':', color=cmap(c))
            for id_ in cur_id2:
                ax.plot((cur_start, cur_end), (id_ + .5, id_ + .5), lw=1.5, ls=':', color=cmap(c))                   

        if i != len(assign_path) -2:
            c = 1 - share_pre_score[i+1] / (share_cur_score[i+1] + 0.001)
            ax.plot((cur_end, cur_end), (0, len(assigned_ref_ids)), color=cmap(c),
                    lw=1)
    for i in range(total):
        plt.text(length + text_pad, i + .5,
                 assigned_ref_gene[i], ha='left', va='center')
    plt.ylim([0, total])
    plt.xlabel('Gene orf position (nt)')
 
    plt.tight_layout()
    plt.savefig('%s.png' % figname)
    plt.close()

    

def draw_gene_assign_dotplot(assign_path, rgene_id: dict, gene: str, length: int,
                     figname: str, qseq: str, rseqs: dict, min_rec: int, top=0, ksize=30):


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

    # Draw self-matrix
    fig = plt.figure(figsize=(8, 8), dpi=300)
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
            sloc.extend([l for l in locs if reg[0] <= l < reg[1]])
        for reg in rec_locations:
            rloc.extend([l for l in locs if reg[0] <= l < reg[1]])
        for reg in homo_locations:
            hloc.extend([l for l in locs if reg[0] <= l < reg[1]])
        self_x.extend(sloc * len(sloc))
        self_y.extend([s for s in sloc for t in sloc])
        rec_x.extend(rloc * len(rloc))
        rec_y.extend([s for s in rloc for t in rloc])
        homo_x.extend(hloc * len(hloc))
        homo_y.extend([s for s in hloc for t in hloc])
        all_x.extend(locs * len(locs))
        all_y.extend([u for u in locs for v in locs])
    plt.scatter(all_x, all_y, s=.5, color='grey')
    plt.scatter(homo_x, homo_y, color='cyan', s=1, label='Homology Arm')
    plt.scatter(self_x, self_y, color='blue', s=1, label='Non-recombined')
    plt.scatter(rec_x, rec_y, color='red', s=1, label='Recombined')
    plt.xlabel('{} (nt)'.format(gene))
    plt.ylabel('{} (nt)'.format(gene))
    plt.legend()
    plt.tight_layout()
    plt.savefig('%s.0.png' % figname)
    plt.close()
    # Draw others
    for pid, (rid, _) in enumerate(ref_rec_length, 1):
        
        rname = assigned_ref_gene[rid]
        rseq = rseqs[rname]
        ref_rec = ref_locations[rid]
        ref_homo = []
        if rid in ref_homo_locations:
            ref_homo = ref_homo_locations[rid]
        # Obtain kmers from query
        rec_kdict = {}
        homo_kdict = {}
        for reg in ref_rec:
            for i in range(reg[0], reg[1]):
                kmer = qseq[i: i+ ksize]
                rec_kdict[kmer] = True
        for reg in ref_homo:
            for i in range(reg[0], reg[1]):
                kmer = qseq[i: i+ ksize]
                homo_kdict[kmer] = True
        kdict = {}
        for i in range(len(rseq) - ksize + 1):
            kmer = rseq[i: i+ksize]
            kdict.setdefault(kmer, [])
            kdict[kmer].append(i)
        rec_x = []
        rec_y = []
        homo_x = []
        homo_y = []
        both_x = []
        both_y = []
        all_x = []
        all_y = []
        for kmer, locs in kdict.items():
            tmpx = locs * len(locs)
            tmpy = [l for l in locs for m in locs]
            all_x.extend(tmpx)
            all_y.extend(tmpy)
            if kmer in rec_kdict and kmer in homo_kdict:
                both_x.extend(tmpx)
                both_y.extend(tmpy)
            elif kmer in rec_kdict:
                rec_x.extend(tmpx)
                rec_y.extend(tmpy)
            elif kmer in homo_kdict:
                homo_x.extend(tmpx)
                homo_y.extend(tmpy)      
        fig = plt.figure(figsize=(8, 8), dpi=300)
        plt.xlim(0, len(rseq))
        plt.ylim(0, len(rseq)) 
        plt.scatter(all_x, all_y, color='grey', s=.5)
        plt.scatter(both_x, both_y, color='orange', s=1, label='both')
        plt.scatter(homo_x, homo_y, color='cyan', s=1, label='Homology Arm')
        plt.scatter(rec_x, rec_y, color='red', s=1, label='Recombined')
        plt.xlabel('{} (nt)'.format(rname))
        plt.ylabel('{} (nt)'.format(rname))
        plt.legend()
        plt.tight_layout()
        plt.savefig('{}.{}.png'.format(figname, pid))
        plt.close()
        
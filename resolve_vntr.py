#!/usr/bin/env python3

import bisect
from collections import  defaultdict
from bam_processing import ReadSegment

def read_vntr_file(vntr_file): ## Add check for chromosome names 
    vntr_list = defaultdict(list)
    vntrs = []
    with open(vntr_file)as f:
        for line in f:
            vntrs.append(line.strip().split())
        for tr in vntrs:
            if not vntr_list[tr[0]]:
                vntr_list[tr[0]].append([int(tr[1])])
                vntr_list[tr[0]].append([int(tr[2])])
            else:
                vntr_list[tr[0]][0].append(int(tr[1]))
                vntr_list[tr[0]][1].append(int(tr[2]))
    return vntr_list


def resolve_vntr_ins(seg,vntr_list):
    tr_reg = vntr_list[seg.ref_id]
    if tr_reg:
        strt = bisect.bisect_right(tr_reg[0],seg.ref_end)
        end = bisect.bisect_left(tr_reg[1],seg.ref_end)
        if strt - end == 1:
            return([(seg.ref_id, tr_reg[0][end]), (seg, '', seg.segment_length)])
    
    
def order_bp(seg1, seg2):
    ref_bp1 = seg1.ref_end if seg1.strand == "+" else seg1.ref_start
    sign1 = 1 if seg1.strand == "+" else -1
    ref_bp2 = seg2.ref_start if seg2.strand == "+" else seg2.ref_end
    sign2 = -1 if seg2.strand == "+" else 1
    if sign1 == sign2:
        return False
    bp_order = [(ref_bp1, sign1, seg1),(ref_bp2, sign2, seg2)]
    bp_length = bp_order[0][0] * bp_order[0][1] + bp_order[1][0] * bp_order[1][1]
    return bp_order, bp_length


def resolve_vntr_split(s1, s2, vntr_list):
    tr_reg = vntr_list[s1.ref_id]
    if tr_reg:
        bp_order, bp_length = order_bp(s1, s2)
        strt = bisect.bisect_left(tr_reg[0],bp_order[0][0])
        end = bisect.bisect_left(tr_reg[1],bp_order[1][0])
        if strt - end == 1:
            return([(s1.ref_id, tr_reg[0][end]),(bp_order[0][2], bp_order[1][2], bp_length)])
                
                
    
def resolve_read_vntr(read, vntr_list, min_sv_size):
    seg_in_vntr = defaultdict(list)
    read.sort(key = lambda s:s.read_start)
    ins_segs = []
    split_segs = []
    new_read = []
    split_seg_vntr = []
    s2 = []
    for s in read:
        if s.is_insertion:
            ins_segs.append(s)
        else:
            split_segs.append(s)
    for ins in ins_segs:
        ins_vntr = resolve_vntr_ins(ins, vntr_list)
        if not ins_vntr:
            new_read.append(ins)
        else:
            seg_in_vntr[ins_vntr[0]].append(ins_vntr[1])
    if len(split_segs) < 2:
        new_read.append(split_segs[0])
    else:
        for s1, s2  in zip(split_segs[:-1], split_segs[1:]):
            if not s1.ref_id == s2.ref_id:
                new_read.append(s1)
                continue
            if not s1.strand == s2.strand: ###Check again
                new_read.append(s1)
                continue
            split_seg_vntr = resolve_vntr_split(s1, s2, vntr_list)
            if split_seg_vntr:
                seg_in_vntr[split_seg_vntr[0]].append(split_seg_vntr[1])
            else:
                new_read.append(s1)
        if not split_seg_vntr and s2:
            new_read.append(s2)
    for key, bp_in_vntr in seg_in_vntr.items():
        bp_len = 0
        segments = []
        for bp in bp_in_vntr:
            bp_len += bp[2]
            segments.append(bp[0])
            if not bp[0].is_insertion:
                segments.append(bp[1])
        segments.sort(key = lambda s:s.ref_start)
        s1 = segments[0]
        s2 = segments[-1]
        if bp_len < -1 * min_sv_size:
            s1.ref_end = key[1]
            s1_segment_length = abs(s1.ref_end - s1.ref_start)
            s2.ref_start = key[1] - bp_len
            if s1.strand == '+':
                s1.read_start = s1.read_start
                s1.read_end = s1.read_start + s1_segment_length
                s2.read_start = s1.read_start + s1_segment_length + 1
                s2.read_end = s2.read_end
            if s1.strand == '-':
                s1.read_start = s1.read_end - s1_segment_length
                s1.read_end = s1.read_end
                s2.read_start = s2.read_start
                s2.read_end = s1.read_end - s1_segment_length+1
            s1.segment_length = abs(s1.read_end - s1.read_start)
            s2.segment_length = abs(s2.read_end - s2.read_start)
            new_read.append(s1)
            new_read.append(s2)
        elif bp_len > min_sv_size:
            s1_segment_length = s2.ref_end - s1.ref_start
            dist = key[1] - s1.ref_start
            mm_rate = 0
            ins_start = s1.read_start + dist
            ins_end = ins_start + bp_len
            new_read.append(ReadSegment(ins_start, ins_end, key[1], key[1], s1.read_id,
                                    s1.ref_id, s1.strand, s1.read_length, bp_len, s1.haplotype, s1.mapq, s1.genome_id, mm_rate, True))
            if split_seg_vntr:
                if s1.strand == '+':
                    s1_read_end = s2.read_end
                    s1_read_start = s1.read_start
                else:
                    s1_read_start = s2.read_start
                    s1_read_end = s1.read_end
                new_read.append(ReadSegment(s1_read_start, s1_read_end, s1.ref_start, s2.ref_end, s1.read_id,s1.ref_id, s1.strand, 
                                            s1.read_length, s1_segment_length, s1.haplotype, s1.mapq, s1.genome_id, s1.mismatch_rate, False))
        else:
            s1_segment_length = s2.ref_end - s1.ref_start
            if s1.strand == '+':
                s1_read_end = s2.read_end
                s1_read_start = s1.read_start
            else:
                s1_read_start = s2.read_start
                s1_read_end = s1.read_end
            new_read.append(ReadSegment(s1_read_start, s1_read_end, s1.ref_start, s2.ref_end, s1.read_id, s1.ref_id, s1.strand, 
                                        s1.read_length, s1_segment_length, s1.haplotype, s1.mapq, s1.genome_id, s1.mismatch_rate, False))
    new_read.sort(key = lambda s:s.read_start)
    return new_read

def resolve_vntr(segments_by_read, vntr_file, min_sv_size):
    vntr_list = read_vntr_file(vntr_file)
    resolved_reads = defaultdict(list)
    for read_id, read in segments_by_read.items():
        resolved_reads[read_id] = resolve_read_vntr(read, vntr_list, min_sv_size)
    return resolved_reads
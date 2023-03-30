#!/usr/bin/env python3

import bisect
from collections import  defaultdict
from bam_processing import ReadSegment, add_read_qual
import logging

logger = logging.getLogger()


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
            return([(seg.ref_id, tr_reg[0][end], tr_reg[1][end]), (seg, '', seg.segment_length)])
    
    
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
        strt = bisect.bisect_right(tr_reg[0],min([bp_order[0][0], bp_order[1][0]]))
        end = bisect.bisect_left(tr_reg[1],max([bp_order[0][0], bp_order[1][0]]))
        if strt - end == 1:
            return([(s1.ref_id, tr_reg[0][end], tr_reg[1][end]),(bp_order[0][2], bp_order[1][2], bp_length)])
        
def filter_vntr_only_segments(split_segs, vntr_list):
    OVERLAP_THR = 0.95
    for s1 in split_segs:
        tr_reg = vntr_list[s1.ref_id]
        if tr_reg:
            strt = bisect.bisect_right(tr_reg[0],s1.ref_start)
            end = bisect.bisect_left(tr_reg[1],s1.ref_end)
            if strt - end == 1:
                vntr_len = tr_reg[1][end] - tr_reg[0][end]
                if s1.segment_length > vntr_len * OVERLAP_THR:
                    return True
        
def calc_new_segments(segments, clipped_segs, vntr_strt, vntr_end, bp_len, split_seg_vntr, min_sv_size):
    seg_span=[]
    is_end = False
    seg_len = 0
    for seg in segments:
        if seg.ref_start >= vntr_strt and seg.ref_end <= vntr_end:
            seg_len += seg.segment_length
            if seg.is_end:
                is_end = True
        else:
            seg_span.append(seg)
    mm_rate = 0
    vntr_len = vntr_end - vntr_strt
    if len(seg_span) == 0:
        return clipped_segs
    if len(seg_span) == 1:
        s1 = seg_span[0]
        if abs(bp_len) < min_sv_size:
            if s1.ref_start > vntr_strt:
                s1.ref_start = vntr_strt
            else:
                s1.ref_end = vntr_end
            s1.segment_length = s1.ref_end - s1.ref_start
            return [s1] + clipped_segs 
        if is_end:
            if bp_len < -1 * min_sv_size:
                return [s1] + clipped_segs
            elif bp_len > min_sv_size:
                new_seg = ReadSegment(0, bp_len, vntr_strt, vntr_strt, s1.read_id,
                                        s1.ref_id, s1.strand, s1.read_length,0, s1.haplotype, s1.mapq, s1.genome_id,mm_rate, False)
                new_seg.is_clipped = True
                return [s1, new_seg]
        else:
            if bp_len > min_sv_size:
                new_seg = ReadSegment(0, bp_len, vntr_strt, vntr_strt, s1.read_id,
                                        s1.ref_id, s1.strand, s1.read_length,bp_len, s1.haplotype, s1.mapq, s1.genome_id,mm_rate, True)
                if s1.ref_start > vntr_strt:
                    s1.ref_start = vntr_strt
                else:
                    s1.ref_end = vntr_end
                s1.segment_length = s1.ref_end - s1.ref_start
                clipped_segs += [s1, new_seg]
                return clipped_segs
            else:
                if s1.ref_start > vntr_strt:
                    s1.ref_start -= seg_len
                else:
                    s1.ref_end += seg_len
                    s1.segment_length += seg_len
                return [s1] + clipped_segs 
    else:
        s1 = seg_span[0]
        s2 = seg_span[-1]
        if bp_len < -1 * min_sv_size:
            s1.ref_end = vntr_strt
            s1.segment_length = abs(s1.ref_end - s1.ref_start)
            s2.segment_length = abs(s2.ref_end - s2.ref_start)
            s2.ref_start = vntr_strt - bp_len
            if s1.strand == '+':
                s1.read_end = s1.read_start + s1.segment_length
                s2.read_start = s1.read_end + 1
            if s1.strand == '-':
                s2.read_end= s2.read_start + s2.segment_length
                s1.read_start = s2.read_end + 1
            #s1.segment_length = abs(s1.read_end - s1.read_start)
            #s2.segment_length = abs(s2.read_end - s2.read_start)
            return [s1, s2] + clipped_segs
        elif bp_len > min_sv_size:
            s1.ref_end = s2.ref_end
            s1.segment_length = s1.ref_end - s1.ref_start
            dist = vntr_strt - s1.ref_start
            mm_rate = 0
            ins_start = s1.read_start + dist
            ins_end = ins_start + bp_len
            s2_new = ReadSegment(ins_start, ins_end, vntr_strt, vntr_strt, s1.read_id,
                                    s1.ref_id, s1.strand, s1.read_length, bp_len, s1.haplotype, s1.mapq, s1.genome_id, mm_rate, True)
            if split_seg_vntr:
                if s1.strand == '+':
                    s1.read_end = s2.read_end
                else:
                    s1.read_start = s2.read_start
            return [s1, s2_new] + clipped_segs 
        else:
            s1.ref_end = s2.ref_end
            s1.segment_length = s1.ref_end - s1.ref_start
            if s1.strand == '+':
                s1.read_end = s2.read_end
            else:
                s1.read_start = s2.read_start
            return [s1] + clipped_segs 

def check_spanning(new_read, vntr_loc):
    vntr_strt = vntr_loc[1]
    vntr_end = vntr_loc[2]
    vntr_ref = vntr_loc[0]
    for seg in new_read:
        if seg.is_clipped:
            continue
        if seg.ref_id == vntr_ref and seg.ref_start < vntr_strt and seg.ref_end > vntr_end:
            return True

def resolve_read_vntr(read, vntr_list, min_sv_size):

    seg_in_vntr = defaultdict(list)
    read.sort(key = lambda s:s.read_start)
    ins_segs = []
    split_segs = []
    clipped_segs = []
    new_read = []
    split_seg_vntr = []
    s2 = []
    for s in read:
        if s.is_insertion:
            ins_segs.append(s)
        elif s.is_clipped:
            clipped_segs.append(s)
        else:
            split_segs.append(s)
            
    if filter_vntr_only_segments(split_segs, vntr_list):
        return new_read
    
    for ins in ins_segs:
        ins_vntr = resolve_vntr_ins(ins, vntr_list)
        if not ins_vntr:
            new_read.append(ins)
        else:
            seg_in_vntr[ins_vntr[0]].append(ins_vntr[1])
            
    if len(split_segs) < 2:
        new_read.append(split_segs[0])
        new_read += clipped_segs
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
        segments = []
        bp_len = 0
        
        for bp in bp_in_vntr:
            bp_len += bp[2]
            if not bp[0].is_insertion:
                segments.append(bp[0])
                segments.append(bp[1])
                
        if not segments:
            s1 = bp[0]
            new_read.append(ReadSegment(s1.read_start + 1, s1.read_start + bp_len, key[1], key[1], s1.read_id,
                                        s1.ref_id, s1.strand, s1.read_length,bp_len, s1.haplotype, s1.mapq, s1.genome_id,0, True))
            if not check_spanning(new_read, key):
                new_read[-1].is_clipped = True
                new_read[-1].is_insertion = False
        else:
            segments.sort(key = lambda s:s.ref_start)
            new_segments = calc_new_segments(segments, clipped_segs, key[1], key[2],bp_len, split_seg_vntr, min_sv_size)
            if new_segments:
                new_read += new_segments
                
    new_read.sort(key = lambda s:s.read_start)
    
    return new_read


def resolve_vntr(segments_by_read, vntr_file, min_sv_size):
    vntr_list = read_vntr_file(vntr_file)
    empty_keys = []
    for read_id, read in segments_by_read.items():
        new_read = resolve_read_vntr(read, vntr_list, min_sv_size)
        if new_read:
            segments_by_read[read_id] = new_read
        else:
            empty_keys.append(read_id)

    for key in empty_keys:
        del segments_by_read[key]


def update_segments_by_read(segments_by_read, ref_lengths, thread_pool, args):
    if args.vntr_file:
        resolve_vntr(segments_by_read, args.vntr_file, args.sv_size)
    logger.info("Annotating reads")
    add_read_qual(segments_by_read, ref_lengths, thread_pool, args.min_mapping_quality, args.max_read_error)

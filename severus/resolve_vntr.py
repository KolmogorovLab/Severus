#!/usr/bin/env python3

import bisect
from collections import  defaultdict
import logging
import numpy as np

from severus.bam_processing import add_read_qual, read_vntr_file, _vntr_segment_status, _vntr_indel_size, VNTR_MIN_INDEL

logger = logging.getLogger()

def order_bp(seg1, seg2):
    ref_bp1 = seg1.ref_end if seg1.strand == 1 else seg1.ref_start
    sign1 = 1 if seg1.strand == 1 else -1
    ref_bp2 = seg2.ref_start if seg2.strand == 1 else seg2.ref_end
    sign2 = -1 if seg2.strand == 1 else 1
    if sign1 == sign2:
        return False
    bp_order = [(ref_bp1, sign1, seg1),(ref_bp2, sign2, seg2)]
    bp_length = bp_order[0][0] * bp_order[0][1] + bp_order[1][0] * bp_order[1][1]
    return bp_order, bp_length

def get_vntr_split(s1, s2, vntr_list):
    tr_reg = vntr_list[s1.ref_id]
    if tr_reg:
        bp_order, bp_length = order_bp(s1, s2)
        strt = bisect.bisect_right(tr_reg[2], min([bp_order[0][0], bp_order[1][0]]))
        end = bisect.bisect_left(tr_reg[3], max([bp_order[0][0], bp_order[1][0]]))
        if strt - end == 1:
            return([(s1.ref_id, tr_reg[0][end], tr_reg[1][end], tr_reg[2][end], tr_reg[3][end], tr_reg[5][end]), [bp_order[0][2], bp_order[1][2]], bp_length])

def filter_vntr_only_segments(split_segs, vntr_list):
    OVERLAP_THR = 0.95
    for s1 in split_segs:
        tr_reg = vntr_list[s1.ref_id]
        if tr_reg:
            strt = bisect.bisect_right(tr_reg[2],s1.ref_start)
            end = bisect.bisect_left(tr_reg[3],s1.ref_end)
            if strt - end == 1:
                vntr_len = tr_reg[1][end] - tr_reg[0][end]
                if s1.segment_length > vntr_len * OVERLAP_THR:
                    s1.is_pass = 'vntr_only'    
      
def resolve_read_vntr(read, vntr_list, min_sv_size):
    def make_single_read(segls):
        s1 = min(segls, key=lambda s: s.ref_start)
        s2 = max(segls, key=lambda s: s.ref_end)
        s2.align_len += s1.align_len
        s2.segment_length += s1.segment_length
        s2.ref_start = min(s2.ref_start, s1.ref_start)
        s2.ref_start_ori = min(s2.ref_start_ori, s1.ref_start_ori)
        s2.ref_end = max(s2.ref_end, s1.ref_end)
        s2.ref_end_ori = max(s2.ref_end_ori, s1.ref_end_ori)
        s2.read_start = min(s1.read_start, s2.read_start) 
        s2.read_end = max(s1.read_end, s2.read_end)
        s2.is_primary = True if s1.is_primary or s2.is_primary else None
        s2.ID = ','.join([str(s1.ID.split(',')[0]), str(s1.ref_start), str(s1.ref_end)])
        s1.is_pass = 'vntr_merge'
        return s2
    def vntr_entry_score(vntr_entry):
        _region_id, n_ins, n_del, _vntr_status, _indel_pos = vntr_entry
        return abs(n_ins - n_del), max(n_ins, n_del), n_ins + n_del
    def merge_vntr_entries(cur_entry, new_entry):
        if cur_entry[3] == 2 or new_entry[3] == 2:
            return max([cur_entry, new_entry], key=vntr_entry_score)
        cur_region_id, cur_ins, cur_del, cur_status, cur_indel_pos = cur_entry
        _region_id, n_ins, n_del, vntr_status, indel_pos = new_entry
        best_pos = max([cur_entry, new_entry], key=vntr_entry_score)[4]
        return (cur_region_id, cur_ins + n_ins, cur_del + n_del, max(cur_status, vntr_status), best_pos or cur_indel_pos or indel_pos)
    def update_segment_vntr(segment, region_id, n_ins, n_del, vntr_status, indel_pos):
        if _vntr_indel_size(n_ins, n_del) <= VNTR_MIN_INDEL:
            return
        vntrs = list(segment.vntr or [])
        new_entry = (region_id, n_ins, n_del, vntr_status, indel_pos)
        for i, cur_entry in enumerate(vntrs):
            cur_region_id = cur_entry[0]
            if cur_region_id == region_id:
                vntrs[i] = merge_vntr_entries(cur_entry, new_entry)
                segment.vntr = vntrs
                return
        vntrs.append(new_entry)
        segment.vntr = vntrs
    seg_in_vntr = defaultdict(list)
    read.sort(key = lambda s:(s.align_start, s.read_start))
    split_segs = []
    new_read = []
    split_seg_vntr = []
    len_in_vntr = defaultdict(int)
    for s in read:
        if not s.is_clipped and not s.is_insertion:
            split_segs.append(s)
    filter_vntr_only_segments(split_segs, vntr_list)
    spanning = [s for s in split_segs if not s.is_pass == 'vntr_only']
    if len(split_segs) > 1 and len(spanning) > 1:
        for s1, s2  in zip(split_segs[:-1], split_segs[1:]):
            if not s1.ref_id == s2.ref_id:
                continue
            if not s1.strand == s2.strand:
                continue
            split_seg_vntr = get_vntr_split(s1, s2, vntr_list)
            if split_seg_vntr:
                seg_in_vntr[split_seg_vntr[0]] += split_seg_vntr[1]
                len_in_vntr[split_seg_vntr[0]] += split_seg_vntr[2]
    seg_in_vntr = dict(sorted(seg_in_vntr.items(), key=lambda x: x[0][1]))
    for key, segls in seg_in_vntr.items():
        merged_segment = make_single_read(segls)
        val = len_in_vntr[key]
        if val > 0:
            n_ins, n_del = val, 0
        else:
            n_ins, n_del = 0, abs(val)
        vntr_status = _vntr_segment_status(merged_segment, key[3], key[4])
        if vntr_status is not None:
            update_segment_vntr(merged_segment, key[5], n_ins, n_del, vntr_status, merged_segment.ref_start)
    new_read = [s for s in read if not s.is_pass == 'vntr_merge' and not s.is_pass == 'vntr_only']
    return new_read


def remove_dedup_segments(segments_by_read):
    
    for i, read in enumerate(segments_by_read):
        dedup_segments_ins = []
        dedup_segments = []
        read.sort(key = lambda s:s.read_start)
        for seg in read:
            if seg.is_insertion or seg.is_clipped:
                dedup_segments_ins.append(seg)
            elif not dedup_segments or dedup_segments[-1].read_start != seg.read_start:
                dedup_segments.append(seg)
        segments_by_read[i] = dedup_segments + dedup_segments_ins
        
            
def resolve_vntr(segments_by_read, vntr_file, min_sv_size):
    vntr_list = read_vntr_file(vntr_file)
    for i, read in enumerate(segments_by_read):
        if not read:
            continue
        new_read = resolve_read_vntr(read, vntr_list, min_sv_size)
        segments_by_read[i] = new_read
            


def update_segments_by_read(segments_by_read, mismatch_histograms, bg_mm, ref_lengths,read_qual,read_qual_len, args):
    remove_dedup_segments(segments_by_read)
    bg_mm = float(np.median(bg_mm))
    if args.vntr_file:
        resolve_vntr(segments_by_read, args.vntr_file, args.sv_size)
    logger.info("Annotating reads")
    add_read_qual(segments_by_read, ref_lengths, bg_mm, mismatch_histograms,read_qual,read_qual_len, args)

    

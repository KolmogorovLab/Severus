import pysam
import numpy as np
import bisect
import gzip
from collections import  defaultdict
import logging

logger = logging.getLogger()
COV_WINDOW_MM  = 1000
COV_WINDOW  = 1000
K_MM = 1000
NUM_HAPLOTYPES = 3
VNTR_MIN_INDEL = 10

class ReadSegment(object):
    __slots__ = ("align_start", "read_start", "read_end", "ref_start", "ref_end","ref_start_ori", "ref_end_ori", "read_id", "ref_id",
                 "strand", "read_length",'align_len','segment_length', "haplotype", "mapq", "genome_id",
                 'mismatch_rate', 'is_pass', "is_insertion", "is_clipped", 'error_rate', 'ins_seq', 'is_primary', 'ins_pos', 'bp_pos', "ID", "PS", "vntr", "in_whitereg")
    def __init__(self, align_start, read_start, read_end, ref_start, ref_end,ref_start_ori, ref_end_ori, read_id, ref_id,
                 strand, read_length,align_len,segment_length, haplotype, mapq, genome_id,
                 mismatch_rate, is_insertion, error_rate, is_primary, PS):
        self.align_start = align_start
        self.read_start = read_start
        self.read_end = read_end
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_start_ori = ref_start_ori
        self.ref_end_ori = ref_end_ori
        self.read_id = read_id
        self.ref_id = ref_id
        self.strand = strand
        self.read_length = read_length
        self.align_len = align_len
        self.segment_length = segment_length
        self.haplotype = haplotype
        self.mapq = mapq
        self.genome_id = genome_id
        self.mismatch_rate = mismatch_rate
        self.is_pass = ''
        self.is_insertion = is_insertion
        self.is_clipped = False
        self.error_rate = error_rate
        self.ins_seq = None
        self.is_primary = is_primary
        self.ins_pos = None
        self.bp_pos = None
        self.ID = None
        self.PS = PS
        self.vntr = None
        self.in_whitereg = False
    def __str__(self):
        return "".join(["read_start=", str(self.read_start), " read_end=", str(self.read_end), " ref_start=", str(self.ref_start),
                         " ref_end=", str(self.ref_end), " read_id=", str(self.read_id), " ref_id=", str(self.ref_id), " strand=", str(self.strand),
                         " read_length=", str(self.read_length), " haplotype=", str(self.haplotype),
                         " mapq=", str(self.mapq), "mismatch_rate=", str(self.mismatch_rate), " read_qual=", str(self.is_pass), " genome_id=", str(self.genome_id)])
    def get_pos(self, bp_pos):
        if bp_pos == "right":
            ref_bp = self.ref_end if self.strand == 1 else self.ref_start
            ref_bp_ori = self.ref_end_ori if self.strand == 1 else self.ref_start_ori
            sign = 1 if self.strand == 1 else -1
        elif bp_pos == "left":
            ref_bp = self.ref_start if self.strand == 1 else self.ref_end
            ref_bp_ori = self.ref_start_ori if self.strand == 1 else self.ref_end_ori
            sign = -1 if self.strand == 1 else 1
        return ref_bp, ref_bp_ori, sign

def _vntr_segment_status(segment, region_start, region_end):
    if segment.ref_start <= region_start and segment.ref_end >= region_end:
        return 2
    if segment.ref_start >= region_start and segment.ref_end <= region_end:
        return 0
    left_spanning = segment.ref_start <= region_start < segment.ref_end <= region_end
    right_spanning = region_start <= segment.ref_start < region_end <= segment.ref_end
    if left_spanning or right_spanning:
        return 1
    return None

def _vntr_indel_size(ins_count, del_count):
    return max(ins_count, del_count)

def _set_segments_vntr(read_segments, vntr, vntr_index):
    if not vntr or not vntr_index:
        return
    region_starts, region_ends, region_ids = vntr_index
    region_bounds = {
        region_id: (region_start, region_end)
        for region_start, region_end, region_id in zip(region_starts, region_ends, region_ids)
    }
    for segment in read_segments:
        segment_vntr = []
        for region_id, ins_count, del_count, indel_pos in vntr:
            if _vntr_indel_size(ins_count, del_count) <= VNTR_MIN_INDEL:
                continue
            region_start, region_end = region_bounds[region_id]
            vntr_status = _vntr_segment_status(segment, region_start, region_end)
            if vntr_status is not None:
                segment_vntr.append((region_id, ins_count, del_count, vntr_status, indel_pos))
        if segment_vntr:
            segment.vntr = segment_vntr

def _region_names_for_pos(pos, region_starts, region_ends, region_names):
    names = []
    region_idx = bisect.bisect_right(region_ends, pos)
    while region_idx < len(region_starts) and region_starts[region_idx] <= pos:
        if pos < region_ends[region_idx]:
            names.append(region_names[region_idx])
        region_idx += 1
    return names

def _set_segments_whitereg(read_segments, white_index):
    if not white_index:
        return
    region_starts, region_ends, region_names = white_index
    if not region_starts:
        return
    for segment in read_segments:
        names = []
        for pos in (segment.ref_start, segment.ref_end):
            names.extend(_region_names_for_pos(pos, region_starts, region_ends, region_names))
        if names:
            segment.in_whitereg = list(dict.fromkeys(names))

def _has_large_vntr_indel(vntr, sv_size):
    return any(max(ins_count, del_count) >= sv_size for _region_id, ins_count, del_count, _indel_pos in vntr)

def _add_full_read_vntr_segment(read_segments, read, align_start, read_length, total_segment_length,
                                haplotype, genome_id, mm_rate, error_rate, is_primary, PS, ID, vntr, vntr_index):
    read_start = align_start
    read_end = min(read_length, read_start + total_segment_length)
    if read.is_reverse:
        read_start, read_end = read_length - read_end, read_length - read_start
    region_starts, region_ends, region_ids = vntr_index
    region_bounds = {
        region_id: (region_start, region_end)
        for region_start, region_end, region_id in zip(region_starts, region_ends, region_ids)
    }
    read_segments.append(ReadSegment(align_start, read_start, read_end, read.reference_start, read.reference_end,
                                     read.reference_start, read.reference_end, read.query_name, read.reference_name,
                                     -1 if read.is_reverse else 1, read_length, total_segment_length, total_segment_length,
                                     haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, is_primary, PS))
    read_segments[-1].ID = ID
    read_segments[-1].vntr = [
        (region_id, ins_count, del_count, vntr_status, indel_pos)
        for region_id, ins_count, del_count, indel_pos in vntr
        if region_id in region_bounds
        if _vntr_indel_size(ins_count, del_count) > VNTR_MIN_INDEL
        for vntr_status in [_vntr_segment_status(read_segments[-1], *region_bounds[region_id])]
        if vntr_status is not None
    ]

def get_segment(read, genome_id,sv_size,use_supplementary_tag, ref_ind, ignore_hp, vntr_regs, control_g, count_window=None, vntr_index=None, white_index=None):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """

    CIGAR_MATCH = [0, 7, 8]
    #CIGAR_MM = 8
    CIGAR_DEL = 2
    CIGAR_INS = 1
    CIGAR_CLIP = [4, 5]
    vntr_sv_size = 10 if control_g else sv_size

    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    MIN_CLIPPED_LENGTH = 500

    ref_start = read.reference_start
    read_segments =[]
    cigar = read.cigartuples
    ID = ','.join([str(ref_ind), str(read.reference_start), str(read.reference_end)])
    vntr_indels = []
    is_primary = None
    if not read.is_supplementary:
        is_primary = True

    num_of_mismatch = 0
    nm = read.get_tag('NM')
    indel = [b for a, b in cigar if a in [CIGAR_INS, CIGAR_DEL]]
    num_of_mismatch = nm - sum(indel)
    total_segment_length = sum([b for a, b in cigar if a not in CIGAR_CLIP + [CIGAR_DEL]])
    if total_segment_length == 0:
        return [[], []] 
    mm_rate = int(num_of_mismatch * K_MM / total_segment_length)
    error_rate = int(nm * K_MM/ total_segment_length)
    read_length = sum([k[1] for k in cigar if k[0] != CIGAR_DEL])
    strand = -1 if read.is_reverse else 1
    sequence = read.query_sequence
    hc = 0
    clp = cigar[-1] if read.is_reverse else cigar[0]
    align_start = clp[1] if clp[0] in CIGAR_CLIP else 0
    
    use_tag = False
    if use_supplementary_tag or is_primary:
        use_tag = True
        
    haplotype = 0
    PS = 0
    if not ignore_hp and use_tag and read.has_tag('HP'):
        haplotype = int(read.get_tag('HP'))
    if not ignore_hp and use_tag and read.has_tag('PS'):
        PS = read.get_tag('PS')
    
    read_inf = []
    vntr_indels = _count_vntr_indels(read, vntr_index, sv_size, count_window)
    
    if (not indel or max(indel) < sv_size) and is_primary and not read.has_tag('SA'):
        if _has_large_vntr_indel(vntr_indels, vntr_sv_size):
            _add_full_read_vntr_segment(read_segments, read, align_start, read_length, total_segment_length,
                                        haplotype, genome_id, mm_rate, error_rate, is_primary, PS, ID, vntr_indels, vntr_index)
        else:
            read_inf = np.array([ref_ind, read.reference_start, read.reference_end,read_length,total_segment_length, haplotype,mm_rate, error_rate, read.mapping_quality], dtype = int)
            
        if cigar[0][0] in CIGAR_CLIP and cigar[0][1] > MIN_CLIPPED_LENGTH:
            pos = read.reference_start if strand == 1 else read.reference_end
            st = -1 if strand == 1 else 1
            read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                        read.reference_name, st, read_length, total_segment_length, cigar[0][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None, PS))
            read_segments[-1].is_clipped = True
        if cigar[-1][0] in CIGAR_CLIP and cigar[-1][1] > MIN_CLIPPED_LENGTH:
            pos = read.reference_start if strand == -1 else read.reference_end
            st = 1 if strand == 1 else -1
            read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                        read.reference_name, st, read_length, total_segment_length, cigar[-1][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None,PS))
            read_segments[-1].is_clipped = True
        _set_segments_whitereg(read_segments, white_index)
        return [read_segments, read_inf]

    for token in cigar:
        op = token[0]
        op_len = token[1]
        if op in CIGAR_CLIP:
            if first_clip:
                read_start = op_len
            hc = op_len if op == 5 else 0
        first_clip = False
        if op in CIGAR_MATCH:
            read_aligned += op_len
            ref_aligned += op_len
        if op == CIGAR_DEL:
            if op_len < sv_size:
                ref_aligned += op_len
            else:
                ref_end = ref_start + ref_aligned
                read_end = read_start + read_aligned
                if read.is_reverse:
                    del_start, del_end = read_length - read_end, read_length - read_start
                else:
                    del_start, del_end = read_start , read_end
                vntr_st = check_vntr(vntr_regs, ref_end) 
                vntr_end = check_vntr(vntr_regs, ref_end+op_len)
                if vntr_st is not None and vntr_end is not None and vntr_st == vntr_end:
                    ref_aligned += op_len
                    continue
                read_segments.append(ReadSegment(align_start, del_start, del_end, ref_start, ref_end, ref_start, ref_end, read.query_name,
                                                 read.reference_name, strand, read_length,total_segment_length,read_aligned,
                                                 haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, is_primary, PS))
                read_segments[-1].ID = ID
                read_start = read_end
                ref_start = ref_end+op_len
                read_aligned = 0
                ref_aligned = 0
        if op == CIGAR_INS:
            if op_len < sv_size:
                read_aligned += op_len
            else:
                ins_start = read_start + read_aligned
                read_aligned += op_len
                ins_pos= ref_start + ref_aligned
                if check_vntr(vntr_regs, ins_pos) is not None:
                    continue
                ins_end = read_start + read_aligned
                read_segments.append(ReadSegment(align_start,ins_start, ins_end, ins_pos, ins_pos, ins_pos, ins_pos, read.query_name,
                                                 read.reference_name, strand, read_length,total_segment_length,op_len, haplotype,
                                                 read.mapping_quality, genome_id, mm_rate, True, error_rate, None,PS))
                ins_seq = sequence[ins_start - hc: ins_end - hc]
                read_segments[-1].ins_seq = ins_seq
                read_segments[-1].ins_pos = (read.reference_start, read.reference_end)
                
    if ref_aligned != 0:
        ref_end = ref_start + ref_aligned
        read_end = read_start + read_aligned
        if read.is_reverse:
            read_start, read_end = read_length - read_end, read_length - read_start
        read_segments.append(ReadSegment(align_start, read_start, read_end, ref_start, ref_end, ref_start, ref_end, read.query_name,
                                         read.reference_name, strand, read_length,total_segment_length,read_aligned, haplotype,
                                         read.mapping_quality, genome_id, mm_rate, False, error_rate, is_primary, PS))
        read_segments[-1].ID = ID
    merge_short_seg(read_segments)
    
        
    if is_primary and not read.has_tag('SA'):
        split = [seg for seg in read_segments if not seg.is_insertion]
        if len(split) == 1 and not _has_large_vntr_indel(vntr_indels, sv_size):
            read_segments = [seg for seg in read_segments if seg.is_insertion]
            read_inf = np.array([ref_ind,  read.reference_start, read.reference_end,read_length,total_segment_length, haplotype,mm_rate, error_rate, read.mapping_quality], dtype = int)
            if cigar[0][0] in CIGAR_CLIP and cigar[0][1] > MIN_CLIPPED_LENGTH:
                pos = read.reference_start if strand == -1 else read.reference_end
                st = -1 if strand == 1 else 1
                read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                            read.reference_name, st, read_length, total_segment_length, cigar[0][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None, PS))
                read_segments[-1].is_clipped = True
            if cigar[-1][0] in CIGAR_CLIP and cigar[-1][1] > MIN_CLIPPED_LENGTH:
                pos = read.reference_start if strand == -1 else read.reference_end
                st = 1 if strand == 1 else -1
                read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                            read.reference_name, st, read_length, total_segment_length, cigar[-1][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None, PS))
                read_segments[-1].is_clipped = True
            
    
    _set_segments_vntr(read_segments, vntr_indels, vntr_index)
    _set_segments_whitereg(read_segments, white_index)
    return [read_segments, read_inf]

def merge_short_seg(read):
    seg_to_remove = []
    DEL_SEG_THR = 500
    INS_THR = 1.5
    INS_DIST_THR = 2000
    read.sort(key=lambda s: s.ref_start)
    init = True
    strand = read[0].strand

    for seg in read:
        if seg.is_insertion:
            continue
        if init or seg.segment_length > DEL_SEG_THR:
            init = False
            init_seg = seg
        elif seg.segment_length <= DEL_SEG_THR:
            init_seg.ref_end = init_seg.ref_end + seg.segment_length
            seg_to_remove.append(seg)

    init = True
    for seg in read:
        if not seg.is_insertion:
            continue
        if init:
            init = False
            init_seg = seg
        elif (seg.ref_start - init_seg.ref_start) > min([INS_THR * (init_seg.segment_length + seg.segment_length), INS_DIST_THR]):
            init_seg = seg
        else:
            init_seg.ref_end = int(np.mean([init_seg.ref_end, seg.ref_end]))
            init_seg.ref_start = init_seg.ref_end
            init_seg.ref_end_ori = init_seg.ref_end
            init_seg.ref_start_ori = init_seg.ref_start
            init_seg.segment_length = init_seg.segment_length + seg.segment_length
            if strand == -1:
                init_seg.read_start = init_seg.read_end - init_seg.segment_length
                init_seg.ins_seq = seg.ins_seq + init_seg.ins_seq
            else:
                init_seg.read_end = init_seg.read_start + init_seg.segment_length
                init_seg.ins_seq += seg.ins_seq
            seg_to_remove.append(seg)

    if seg_to_remove:
        for seg in seg_to_remove:
            read.remove(seg)


def extract_clipped_end(segments_by_read):
    MIN_CLIPPED_LENGTH = 500
    for read in segments_by_read:
        read2 = [seg for seg  in read if not seg.is_insertion and seg.is_pass == 'PASS']
        if not read2:
            continue
        read2.sort(key=lambda s: s.read_start)
        s1 = read2[0]
        s2 = read2[-1]
        if s1.read_start > MIN_CLIPPED_LENGTH:
            pos = s1.ref_start if s1.strand == 1 else s1.ref_end
            st = -1 if s1.strand == 1 else 1
            read.append(ReadSegment(0, 0, s1.read_start, pos, pos, pos, pos, s1.read_id,
                                    s1.ref_id, st, s1.read_length, s1.align_len, s1.segment_length, s1.haplotype, s1.mapq, s1.genome_id, s1.mismatch_rate, False, s1.error_rate, None, s1.PS))
            read[-1].is_clipped = True
            read[-1].is_pass = 'PASS'
        end_clip_length = s2.read_length - s2.read_end
        if end_clip_length > MIN_CLIPPED_LENGTH:
            pos = s2.ref_end if s2.strand == 1 else s2.ref_start
            st = 1 if s2.strand == 1 else -1
            read.append(ReadSegment(s2.read_end, s2.read_end, s2.read_length, pos, pos, pos, pos, s2.read_id,
                                    s2.ref_id, st, s2.read_length, s2.align_len, s2.segment_length, s2.haplotype, s2.mapq, s2.genome_id, s2.mismatch_rate, False, s2.error_rate, None, s2.PS))
            read[-1].is_clipped = True
            read[-1].is_pass = 'PASS'
        read.sort(key=lambda s: s.read_start)
        

def get_cov(bam_file, genome_id, ref_id, poslist, min_mapq, ignore_hp, ref = None):

    region_start, region_end = max(0, poslist[0]-3), poslist[-1]+3
    if bam_file.endswith(".bam"):
        aln_file = pysam.AlignmentFile(bam_file, "rb")
            
    elif bam_file.endswith(".cram"):
        aln_file = pysam.AlignmentFile(bam_file, "rc", reference_filename=ref)
    
    cov_list = defaultdict(list)
    BUFF = 50
    BUFF2 = 5
    for pos in poslist:
        cov_list[(ref_id, pos)] = [0,0,0,0,0,0]
    
    ref_lengths = dict(zip(aln_file.references, aln_file.lengths))
        
    if not ref_id in ref_lengths.keys():
        return cov_list
    
    #region_end = min(region_end, ref_lengths[ref_id])
        
    for aln in aln_file.fetch(ref_id, region_start, region_end,  multiple_iterators=True):
        if not aln.is_unmapped:
            if aln.is_secondary:
                strt = bisect.bisect_left(poslist, aln.reference_start)
                end = bisect.bisect_right(poslist, aln.reference_end)
                if not strt == end:
                    for pos in poslist[strt:end]:
                        cov_list[(ref_id, pos)][5]+=1
                        if abs(aln.reference_start - pos) <= BUFF:
                            cov_list[(ref_id, pos)][4]+=1
                            
                        elif abs(aln.reference_end - pos) <= BUFF:
                            cov_list[(ref_id, pos)][3]+=1
                            
            elif aln.mapping_quality > min_mapq:
                strt = bisect.bisect_right(poslist, aln.reference_start)
                end = bisect.bisect_left(poslist, aln.reference_end)
                if not strt == end:
                    if ignore_hp:
                        hp = 0
                    elif aln.has_tag('HP'):
                        hp = int(aln.get_tag('HP'))
                    else:
                        hp = 0
                    for pos in poslist[strt:end]:
                        cov_list[(ref_id, pos)][hp]+=1
                
    return cov_list



def get_coverage_parallel(bam_files, genome_ids, thread_pool, min_mapq, double_breaks, ignore_hp, ref):

    db_list = defaultdict(list)
    covlist = defaultdict(list)
    for db in double_breaks:
        db_list[db.bp_1.ref_id].append(db.bp_1.position)
        db_list[db.bp_2.ref_id].append(db.bp_2.position)
    CHUNK_SIZE = 10000000
    for key, lst in db_list.items():
        lst = list(set(lst))
        lst.sort()
        indices = [i + 1 for (x, y, i) in zip(lst, lst[1:], range(len(lst))) if CHUNK_SIZE < abs(x - y)]
        db_list[key] = [lst[start:end] for start, end in zip([0] + indices, indices + [len(lst)])]
    for genome_id in genome_ids:
        covlist = defaultdict(list)
        tasks = [(bam_files[genome_id], genome_id, ref_id, pos, min_mapq, ignore_hp, ref) for ref_id, poslist in db_list.items() for pos in poslist]

        parsing_results = None
        parsing_results = thread_pool.starmap(get_cov, tasks)
        for item in parsing_results:
            for key, value in item.items():
                covlist[key] = value
        for db in double_breaks:
            db.bp_1.spanning_reads[genome_id] = covlist[(db.bp_1.ref_id, db.bp_1.position)]
            db.bp_2.spanning_reads[genome_id] = covlist[(db.bp_2.ref_id, db.bp_2.position)]
    

def _bed_overlaps(start, end, region_starts, region_ends, region_ids):
    if end <= start:
        return []
    overlaps = []
    region_idx = bisect.bisect_right(region_ends, start)
    while region_idx < len(region_starts) and region_starts[region_idx] < end:
        region_start = region_starts[region_idx]
        region_end = region_ends[region_idx]
        overlap = max(0, min(end, region_end) - max(start, region_start))
        if overlap:
            overlaps.append((region_ids[region_idx], overlap, max(start, region_start), region_start, region_end))
        region_idx += 1
    return overlaps

def _has_bed_overlap(start, end, region_starts, region_ends):
    if end <= start:
        return False
    region_idx = bisect.bisect_right(region_ends, start)
    return region_idx < len(region_starts) and region_starts[region_idx] < end

def _region_info_for_pos(pos, region_starts, region_ends, region_ids):
    region_idx = bisect.bisect_right(region_starts, pos) - 1
    if region_idx >= 0 and pos < region_ends[region_idx]:
        return region_ids[region_idx], region_starts[region_idx], region_ends[region_idx]
    return None

def _index_vntr_regions(vntr_regs):
    if not vntr_regs:
        return [], [], []
    return vntr_regs[2], vntr_regs[3], vntr_regs[5]

def _index_region_names(regions):
    if not regions:
        return [], [], []
    return regions[0], regions[1], regions[4]

def _update_vntr_indel_pos(read_indel, indel_len, indel_pos, sv_size):
    if indel_len > sv_size and read_indel[4] is None:
        read_indel[4] = indel_pos
    if indel_len > read_indel[2]:
        read_indel[2] = indel_len
        read_indel[3] = indel_pos

def _count_vntr_indels(read, vntr_index, sv_size, count_window=None):
    if vntr_index is None or read.reference_end is None:
        return []

    region_starts, region_ends, region_ids = vntr_index
    if not region_starts:
        return []

    window_start, window_end = count_window if count_window else (read.reference_start, read.reference_end)
    read_start = max(read.reference_start, window_start)
    read_end = min(read.reference_end, window_end)
    if not _has_bed_overlap(read_start, read_end, region_starts, region_ends):
        return []

    cigar_ins = 1
    cigar_del = 2
    ref_consuming_ops = {0, 2, 3, 7, 8}
    LEN_THR = 3
    read_indels = {}
    ref_pos = read.reference_start

    for op, op_len in read.cigartuples or []:
        if op == cigar_ins and op_len >= LEN_THR:
            region_info = _region_info_for_pos(ref_pos, region_starts, region_ends, region_ids)
            if window_start <= ref_pos < window_end and region_info is not None:
                region_id, _region_start, _region_end = region_info
                read_indel = read_indels.setdefault(region_id, [0, 0, 0, None, None])
                read_indel[0] += op_len
                _update_vntr_indel_pos(read_indel, op_len, ref_pos, sv_size)
        elif op == cigar_del and op_len >= LEN_THR:
            del_start = ref_pos
            del_end = ref_pos + op_len
            for region_id, overlap, overlap_start, region_start, region_end in _bed_overlaps(max(del_start, window_start), min(del_end, window_end), region_starts, region_ends, region_ids):
                if del_start <= region_start and del_end >= region_end:
                    continue
                read_indel = read_indels.setdefault(region_id, [0, 0, 0, None, None])
                read_indel[1] += overlap
                _update_vntr_indel_pos(read_indel, overlap, overlap_start, sv_size)

        if op in ref_consuming_ops:
            ref_pos += op_len

    return [
        (region_id, indel[0], indel[1], indel[4] if indel[4] is not None else indel[3])
        for region_id, indel in read_indels.items()
    ]


def get_all_reads(bam_file, region, genome_id,sv_size,use_supplementary_tag, ignore_hp, control_g, ref=None, vntrs=None, white_regs=None):

    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """

    alignments = []
    read_info_final = []
    ncol= 10000
    read_info = np.zeros((ncol,9), dtype = int)
    ref_ind, ref_id, region_start, region_end = region
    vntr_index = _index_vntr_regions(vntrs)
    white_index = _index_region_names(white_regs)
    if bam_file.endswith(".bam"):
        aln_file = pysam.AlignmentFile(bam_file, "rb")
            
    elif bam_file.endswith(".cram"):
        aln_file = pysam.AlignmentFile(bam_file, "rc", reference_filename=ref)
        
    t=0
    for aln in aln_file.fetch(ref_id, region_start, region_end,  multiple_iterators=True):
        if not aln.is_secondary and not aln.is_unmapped:
            if aln.reference_end is None:
                continue
            new_segment, read_inf = get_segment(aln, genome_id, sv_size, use_supplementary_tag, ref_ind, ignore_hp, vntrs, control_g, (region_start, region_end), vntr_index, white_index)
            if new_segment:
                alignments += new_segment
            if not len(read_inf):
                continue
            read_info[t] = read_inf
            t+=1
            if t == ncol:
                t = 0
                if not len(read_info_final):
                    read_info_final = read_info
                else:
                    read_info_final = np.concatenate((read_info_final, read_info), axis = 0)
                read_info = np.zeros((ncol,9), dtype = int)
               
    if not len(read_info_final):
        t = 1 if t ==0 else t
        read_info_final = read_info[0:t-1]
    else:
        read_info_final = np.concatenate((read_info_final, read_info[0:t-1]), axis=0)                
    return (alignments, read_info_final)

def read_whitelist_file(bed_file):
    """Read and merge whitelist BED intervals in the region-reader format."""
    intervals = defaultdict(list)
    if not bed_file:
        return defaultdict(list)
    opener = gzip.open if str(bed_file).endswith('.gz') else open
    with opener(bed_file, 'rt') as bed:
        for line_number, line in enumerate(bed, 1):
            fields = line.split()
            if not fields or fields[0].startswith('#') or fields[0] in ('track', 'browser'):
                continue
            if len(fields) < 3:
                raise ValueError(f'Whitelist BED line {line_number} requires three columns')
            chrom, start, end = fields[:3]
            start, end = int(start), int(end)
            if start < 0 or end < start:
                raise ValueError(f'Invalid whitelist interval on line {line_number}')
            intervals[chrom].append((start, end))
    regions = defaultdict(list)
    region_id = 0
    for chrom, entries in intervals.items():
        merged = []
        for start, end in sorted(entries):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        starts = [x[0] for x in merged]
        ends = [x[1] for x in merged]
        regions[chrom] = [starts, ends, starts[:], ends[:],
                          [f'{chrom}_{start}' for start in starts],
                          list(range(region_id, region_id + len(merged)))]
        region_id += len(merged)
    return regions


def read_vntr_file(vntr_file, bp_tol=50):
    vntr_list = defaultdict(list)
    region_idx = 0
    if not vntr_file:
        return vntr_list
    if vntr_file.endswith(".bed.gz"):
        f = gzip.open(vntr_file, "rt")
    elif vntr_file.endswith(".bed"):
        f = open(vntr_file)
    else:
        raise ValueError("VNTR file must end with .bed or .bed.gz")

    with f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            tr = line.strip().split()
            chrom = tr[0]
            start = int(tr[1])
            end = int(tr[2])
            if end <= start:
                continue
            if not vntr_list[chrom]:
                vntr_list[chrom] = [[], [], [], [], [], []]
            vntr_list[chrom][0].append(start)
            vntr_list[chrom][1].append(end)
            vntr_list[chrom][2].append(start - bp_tol)
            vntr_list[chrom][3].append(end + bp_tol)
            vntr_list[chrom][4].append(f"{chrom}_{start}")
            vntr_list[chrom][5].append(region_idx)
            region_idx += 1
    return vntr_list

def get_vntrs_in_region(vntr_list, chrom, start, end):
    if not vntr_list:
        return []
    if chrom not in vntr_list:
        return []
    tr_reg = vntr_list[chrom] if isinstance(vntr_list, defaultdict) else vntr_list
    if not tr_reg:
        return []

    start_idx = bisect.bisect_left(tr_reg[3], start)
    end_idx = bisect.bisect_right(tr_reg[2], end)
    if start_idx >= end_idx:
        return []
    return [values[start_idx:end_idx] for values in tr_reg]

def check_vntr(vntr_regs, pos):
    if not vntr_regs:
        return None
    start_idx = bisect.bisect_right(vntr_regs[2], pos)
    end_idx = bisect.bisect_left(vntr_regs[3], pos)
    if start_idx - end_idx == 1:
        return end_idx
    return None

def get_all_reads_parallel(bam_file, thread_pool, ref_lengths, genome_id, control_g,
                           coverage_histograms, mismatch_histograms, n90ls, bg_mmls,read_qual,read_qual_len, args):

    CHUNK_SIZE = 10000000
    sv_size = args.sv_size
    use_supplementary_tag = args.use_supplementary_tag
    ignore_hp = args.ignore_hp
    vntr_list = read_vntr_file(getattr(args, 'vntr_file', None))
    white_reg_file = getattr(args, 'white_reg', None) or getattr(args, 'whitelist', None)
    white_reg_list = getattr(args, 'whitelist_regions', None)
    if white_reg_list is None:
        white_reg_list = read_whitelist_file(white_reg_file)
    
    ref = args.ref
    if bam_file.endswith(".bam"):
        all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
            
    elif bam_file.endswith(".cram"):
        all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rc", reference_filename=args.ref).references]

    
    fetch_list = []
    for j, ctg in enumerate(ref_lengths):
        if ctg not in all_reference_ids:
            continue
        ctg_len = ref_lengths[ctg]
        for i in range(0, max(ctg_len // CHUNK_SIZE, 1)):
            reg_start = i * CHUNK_SIZE
            reg_end = (i + 1) * CHUNK_SIZE
            if ctg_len - reg_end < CHUNK_SIZE:
                reg_end = ctg_len
            fetch_list.append((j, ctg, reg_start, reg_end))

    tasks = [
        (bam_file, region, genome_id, sv_size, use_supplementary_tag, ignore_hp, control_g, ref,
         get_vntrs_in_region(vntr_list, region[1], region[2], region[3]),
         get_vntrs_in_region(white_reg_list, region[1], region[2], region[3]))
        for region in fetch_list
    ]
    parsing_results = None
    parsing_results = thread_pool.starmap(get_all_reads, tasks)
    segments_by_read = defaultdict(list)
    for alignments in parsing_results:
        for aln in alignments[0]:
            if aln:
                segments_by_read[aln.read_id].append(aln)
    n90, bg_mm = calc_read_qual(parsing_results, segments_by_read, mismatch_histograms, coverage_histograms, genome_id, ref_lengths, read_qual, read_qual_len, args)
    n90ls.append(n90)
    bg_mmls.append(bg_mm)
    
    return list(segments_by_read.values())


def calc_read_qual(parsing_results,segments_by_read, mismatch_histograms, coverage_histograms, genome_id, ref_lengths, read_qual, read_qual_len, args):
    bg_mm = background_mm_rat(parsing_results, args.multisample)
    update_mm_hist(parsing_results, mismatch_histograms, ref_lengths)
    n90 = get_read_statistics(parsing_results, segments_by_read)
    n90 = min(n90, args.min_aligned_length) if not args.multisample else args.min_aligned_length
    update_cov_hist(parsing_results, coverage_histograms, genome_id, ref_lengths, bg_mm, n90, read_qual, read_qual_len, args)
    return n90, bg_mm

def init_hist(genome_ids, ref_lengths):
    coverage_histograms = {}
    for genome_id in genome_ids:
        for chr_id, chr_len in ref_lengths.items():
            for hp in range(0, NUM_HAPLOTYPES):
                coverage_histograms[(genome_id, hp, chr_id)] = [0 for _ in range(chr_len // COV_WINDOW + 1)]
    return coverage_histograms


def update_cov_hist(parsing_results, coverage_histograms, genome_id, ref_lengths, bg_mm, n90, read_qual, read_qual_len, args):
    for alignments in parsing_results:
        if alignments[1].size == 0:
            continue
        chr_id = list(ref_lengths.keys())[int(alignments[1][0][0])]
        for j in range(len(alignments[1])):
            if alignments[1][j][8] > args.min_mapping_quality and alignments[1][j][6] < bg_mm and alignments[1][j][3] > n90 and abs(alignments[1][j][3] - alignments[1][j][4]) < alignments[1][j][4]:
                read_qual['PASS'] += 1
                read_qual_len['PASS'] += alignments[1][j][4]
                hist_start = alignments[1][j][1] // COV_WINDOW
                hist_end = alignments[1][j][2]// COV_WINDOW
                for i in range(hist_start, hist_end + 1):
                    coverage_histograms[(genome_id, alignments[1][j][5], chr_id)][i] += 1
            else:
                if alignments[1][j][8] < args.min_mapping_quality:
                    read_qual['FAIL_MAPQ'] += 1
                    read_qual_len['FAIL_MAPQ'] += alignments[1][j][4]
                elif alignments[1][j][6] > bg_mm:
                    read_qual['FAIL_MM'] += 1
                    read_qual_len['FAIL_MM'] += alignments[1][j][4]
                elif alignments[1][j][3] < n90:
                    read_qual['FAIL_READLEN'] += 1
                    read_qual_len['FAIL_READLEN'] += alignments[1][j][4]
                elif abs(alignments[1][j][3] - alignments[1][j][4]) > alignments[1][j][4]:
                    read_qual['FAIL_ALNLEN'] += 1
                    read_qual_len['FAIL_ALNLEN'] += alignments[1][j][4]

def update_coverage_hist(coverage_histograms,genome_ids, ref_lengths, segments_by_read, control_genomes, target_genomes):
    
    for read in segments_by_read:
        for seg in read:
            if seg.is_pass == 'PASS' and not seg.is_insertion and not seg.is_clipped:
                hist_start = seg.ref_start_ori // COV_WINDOW
                hist_end = min([seg.ref_end_ori, ref_lengths[seg.ref_id]])// COV_WINDOW
                for i in range(hist_start + 1, hist_end):
                    coverage_histograms[(seg.genome_id, seg.haplotype, seg.ref_id)][i] += 1
    median_cov = []               
    for genome_id in genome_ids:
        by_hp = {}
        for hp in range(0, NUM_HAPLOTYPES):
            by_hp[hp] = []
            for chr_id, _chr_len in ref_lengths.items():
                by_hp[hp].extend(coverage_histograms[(genome_id, hp, chr_id)])

        hp1_cov, hp2_cov, hp0_cov = np.median(by_hp[1]), np.median(by_hp[2]), np.median(by_hp[0])
        logger.info(f"\tMedian coverage by PASS reads for {genome_id} (H1 / H2 / H0): {hp1_cov} / {hp2_cov} / {hp0_cov}")
        median_cov.append(int(max([hp1_cov, hp2_cov, hp0_cov])))

    return max(median_cov)
        

def add_read_qual(segments_by_read, ref_lengths, bg_mm, mismatch_histograms,read_qual,read_qual_len, args):
    min_mapq = args.min_mapping_quality
    min_aligned_length = args.min_aligned_length
    write_segdups_out = args.write_segdups_out
    white_reg_file = getattr(args, 'white_reg', None) or getattr(args, 'whitelist', None)
    white_reg = getattr(args, 'whitelist_regions', None)
    if white_reg is None:
        white_reg = read_whitelist_file(white_reg_file)
    
    mm_hist_high = background_mm_hist(segments_by_read, mismatch_histograms, bg_mm, ref_lengths)

    if write_segdups_out:
        extract_segdups(mm_hist_high, write_segdups_out) ###
    for alignments in segments_by_read:
        if not alignments:
            continue
        label_reads(alignments, min_mapq, bg_mm, mm_hist_high, min_aligned_length, args.multisample, args.min_len, white_reg)

    write_readqual(segments_by_read, args.outpath_readqual,read_qual,read_qual_len)

#Salute to AmberVerhasselt

def in_whitelist(white_reg, ref_id, start, end):
    if not white_reg:
        return None
    if ref_id not in white_reg:
        return None
    regions = white_reg[ref_id]
    starts, ends = regions[0], regions[1]
    return bisect.bisect_right(starts, end) > bisect.bisect_left(ends, start)


def label_reads(read, min_mapq, bg_mm, mm_hist_high, min_aligned_length, multisample, min_len, white_reg, enforce_mapq=False):
    MIN_ALIGNED_RATE = 0 if not min_len == -1 else 0.5
    MIN_ALIGNED_LEN = min_len if not min_len == -1 else min_aligned_length
    #MIN_ALIGNED_LEN = min_aligned_length if multisample else 2000

    for seg in read:
        whitelisted = in_whitelist(white_reg, seg.ref_id, seg.ref_start_ori, seg.ref_end_ori)
        if whitelisted:
            seg.is_pass = 'PASS'
            continue
        if seg.mapq < min_mapq and not enforce_mapq:
            seg.is_pass += '_LOW_MAPQ'
        if high_mm_check(mm_hist_high, bg_mm, seg):
            seg.is_pass += '_HIGH_MM_rate'

    seg_ins = [1 for seg in read if not seg.is_insertion and not seg.is_clipped]
    aligned_len = sum([seg.read_end - seg.read_start for seg in read if not seg.is_clipped]) if seg_ins else read[0].align_len
    aligned_ratio = aligned_len/read[0].read_length

    if aligned_ratio < MIN_ALIGNED_RATE or aligned_len < MIN_ALIGNED_LEN:
        for seg in read:
            if seg.is_pass == 'PASS':
                continue
            seg.is_pass += '_LOW_ALIGNED_LEN'#

    for seg in read:
        if not seg.is_pass:
            seg.is_pass = 'PASS'

def write_readqual(segments_by_read, outpath, read_qual, read_qual_len):
    for read in segments_by_read:
        for seg in read:
            if seg.is_clipped or seg.is_insertion or 'vntr_only' in seg.is_pass:
                continue
            read_qual[seg.is_pass] += 1
            read_qual_len[seg.is_pass] += seg.ref_end - seg.ref_start

    f = open(outpath, "w")
    f.write('Number of segments:')
    f.writelines('{}\t{}\n'.format(k,v) for k, v in read_qual.items())
    f.write('Total length of segments:')
    f.writelines('\n{}\t{}'.format(k,v) for k, v in read_qual_len.items())



def background_mm_rat(parsing_results, multisample):
    COV_WINDOW_BG_MM = 2000
    QT = 0.95 if not multisample else 0.975
    mm_list = []
    for alignments in parsing_results:
        for aln in alignments[1]:
            n_mm = (aln[3] // COV_WINDOW_BG_MM) +1
            mm_list += [aln[6]] * n_mm
        for seg in alignments[0]:
            if seg.is_clipped or seg.is_insertion:
                continue
            n_mm = ((seg.ref_end - seg.ref_start) // COV_WINDOW_BG_MM) +1
            mm_list += [seg.mismatch_rate] * n_mm
    bg_mm = np.quantile(mm_list, QT)
    return bg_mm

def init_mm_hist(ref_lengths):
    mismatch_histograms = defaultdict(list)
    for chr_id, chr_len in ref_lengths.items():
        mismatch_histograms[chr_id] = [[] for _ in range(chr_len // COV_WINDOW_MM + 2)]
    return mismatch_histograms

def update_mm_hist(parsing_results, mismatch_histograms, ref_lengths):
    ref_ind = list(ref_lengths.keys())
    for alignments in parsing_results:
        if alignments[1].size == 0:
            continue
        chr_id = ref_ind[int(alignments[1][0][0])]
        for j in range(len(alignments[1])):
            hist_start = alignments[1][j][1] // COV_WINDOW_MM
            hist_end = alignments[1][j][2] // COV_WINDOW_MM
            for i in range(hist_start, hist_end + 1):
                mismatch_histograms[chr_id][i].append(int(alignments[1][j][6]))
    
def background_mm_hist(segments_by_read, mismatch_histograms, bg_mm, ref_lengths):
    MED_PER = 0.1
    MIN_READ = 5
    MAX_COV = 15
    
    mm_hist_high = {}

    for chr_id, chr_len in ref_lengths.items():
        mm_hist_high[chr_id] = [0 for _ in range(chr_len // COV_WINDOW_MM + 2)]
    
    for read in segments_by_read:
        for seg in read:
            if seg.is_clipped or seg.is_insertion:
                continue
            hist_start = seg.ref_start_ori // COV_WINDOW_MM
            hist_end = min([seg.ref_end_ori, ref_lengths[seg.ref_id]])// COV_WINDOW_MM
            for i in range(hist_start, hist_end + 1):
                mismatch_histograms[seg.ref_id][i].append(seg.mismatch_rate)
                
    for chr_id, mm_rate in mismatch_histograms.items():
        for i, mm_list in enumerate(mm_rate):
            if not mm_list or len(mm_list) < MIN_READ:
                mm_hist_high[chr_id][i] = 0
                continue
            mm_list.sort()
            if mm_list[-4] < bg_mm:
                continue
            med_thr = max([MIN_READ, min([int(len(mm_list)*MED_PER), MAX_COV])])
            if mm_list[med_thr - 1] < bg_mm:
                for k, mm in enumerate(mm_list):
                        if mm > bg_mm:
                            mm_hist_high[chr_id][i] = 1
                            break
    return mm_hist_high
            
 
def extract_segdups(mm_hist_high, write_segdups_out):
    high_mm_region = defaultdict(list)
    for key, mm_high_chrom in mm_hist_high.items():
        for ind, a in enumerate(mm_high_chrom):
            if a == 0:
                continue
            if not high_mm_region[key]:
                high_mm_region[key].append([ind * COV_WINDOW_MM])
                high_mm_region[key].append([(ind + 1) * COV_WINDOW_MM])
            elif ind * COV_WINDOW_MM == high_mm_region[key][1][-1]:
                high_mm_region[key][1][-1]=(ind + 1) * COV_WINDOW_MM
            else:
                high_mm_region[key][0].append(ind * COV_WINDOW_MM)
                high_mm_region[key][1].append((ind + 1) * COV_WINDOW_MM)

    for chr_id, values in high_mm_region.items():
        for (st, end) in zip(values[0], values[1]):
            write_segdups_out.write('\t'.join([chr_id, str(st), str(end)]) + '\n')

def high_mm_check(mm_hist_high, bg_mm, seg):
    strt = seg.ref_start_ori // COV_WINDOW_MM
    end = min([seg.ref_end_ori // COV_WINDOW_MM , len(mm_hist_high[seg.ref_id])-1])
    if any(mm_hist_high[seg.ref_id][strt:end+1]) and seg.mismatch_rate >= bg_mm:
        return True

def _calc_nx(lengths, norm_len, rate):
    n50 = 0
    sum_len = 0
    l50 = 0
    for l in sorted(lengths, reverse=True):
        sum_len += l
        l50 += 1
        if sum_len > rate * norm_len:
            n50 = l
            break
    return l50, n50

def get_read_statistics(parsing_results,segments_by_read):
    ncol= sum([len(a[1]) for a in parsing_results]) + len(segments_by_read)
    read_lengths = np.zeros(ncol, dtype = int)
    alignment_lengths = np.zeros(ncol, dtype = int)
    aln_error = np.zeros(ncol, dtype = int)
    aln_mm = np.zeros(ncol, dtype = int)
    t = 0
    for alignments in parsing_results:
        if alignments[1].size == 0:
            continue
        for i in range(len(alignments[1])):
            alignment_lengths[t] = alignments[1][i][4]
            read_lengths[t] = alignments[1][i][3]
            aln_error[t] = alignments[1][i][7]
            aln_mm[t] = alignments[1][i][6]
            t+=1
            
    for read in segments_by_read.values():
        aln_len = 0
        for aln in read:
            if aln.is_clipped or aln.is_insertion or 'vntr_only' in aln.is_pass:
                continue
            aln_len += aln.segment_length
            
        read_lengths[t] = read[0].read_length
        aln_error[t] = read[0].error_rate
        aln_mm[t] = read[0].mismatch_rate
        alignment_lengths[t] = aln_len
        t+=1
            
    if alignment_lengths.size == 0:
        return None
    
    reads_total_len = np.sum(read_lengths)
    alignment_total_len = np.sum(alignment_lengths)
    aln_rate = alignment_total_len / reads_total_len
    _l50, reads_n50 = _calc_nx(read_lengths, reads_total_len, 0.50)
    _l90, reads_n90 = _calc_nx(read_lengths, reads_total_len, 0.90)
    _l50, aln_n50 = _calc_nx(alignment_lengths, alignment_total_len, 0.50)
    _l50, aln_n90 = _calc_nx(alignment_lengths, alignment_total_len, 0.90)
    [error_25, error_50, error_75] = np.quantile(aln_error, [0.25, 0.50, 0.75])/K_MM
    mm_25, mm_50, mm_75 = np.quantile(aln_mm, [0.25, 0.50, 0.75])/K_MM

    logger.info(f"\tTotal read length: {reads_total_len}")
    logger.info(f"\tTotal aligned length: {alignment_total_len} ({aln_rate:.2f})")
    logger.info(f"\tRead N50 / N90: {reads_n50} / {reads_n90}")
    logger.info(f"\tAlignments N50 / N90: {aln_n50} / {aln_n90}")
    logger.info(f"\tRead error rate (Q25 / Q50 / Q75): {error_25:.4f} / {error_50:.4f} / {error_75:.4f}")
    logger.info(f"\tRead mismatch rate (Q25 / Q50 / Q75): {mm_25:.4f} / {mm_50:.4f} / {mm_75:.4f}")
    return reads_n90

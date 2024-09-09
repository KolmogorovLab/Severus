import pysam
import numpy as np
import bisect
from collections import  defaultdict
import logging
import datetime

logger = logging.getLogger()
COV_WINDOW_MM  = 1000
COV_WINDOW  = 1000
K_MM = 1000
NUM_HAPLOTYPES = 3

class ReadSegment(object):
    __slots__ = ("align_start", "read_start", "read_end", "ref_start", "ref_end","ref_start_ori", "ref_end_ori", "read_id", "ref_id",
                 "strand", "read_length",'align_len','segment_length', "haplotype", "mapq", "genome_id",
                 'mismatch_rate', 'is_pass', "is_insertion", "is_clipped", 'error_rate', 'ins_seq', 'is_primary', 'ins_pos', 'bp_pos')
    def __init__(self, align_start, read_start, read_end, ref_start, ref_end,ref_start_ori, ref_end_ori, read_id, ref_id,
                 strand, read_length,align_len,segment_length, haplotype, mapq, genome_id,
                 mismatch_rate, is_insertion, error_rate, is_primary):
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
    
def get_segment(read, genome_id,sv_size,use_supplementary_tag, ref_ind):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """

    CIGAR_MATCH = [0, 7, 8]
    #CIGAR_MM = 8
    CIGAR_DEL = 2
    CIGAR_INS = 1
    CIGAR_CLIP = [4, 5]

    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    MIN_CLIPPED_LENGTH = 500

    ref_start = read.reference_start
    read_segments =[]
    cigar = read.cigartuples
    
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

    if use_tag and read.has_tag('HP'):
        haplotype = read.get_tag('HP')
    else:
        haplotype = 0
    
    read_inf = []
    
    if (not indel or max(indel) < sv_size) and is_primary and not read.has_tag('SA'):
        read_inf = np.array([ref_ind, read.reference_start, read.reference_end,read_length,total_segment_length, haplotype,mm_rate, error_rate, read.mapping_quality], dtype = int)
        if cigar[0][0] in CIGAR_CLIP and cigar[0][1] > MIN_CLIPPED_LENGTH:
            pos = read.reference_start if strand == 1 else read.reference_end
            read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                        read.reference_name, strand, read_length, total_segment_length, cigar[0][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None))
            read_segments[-1].is_clipped = True
        if cigar[-1][0] in CIGAR_CLIP and cigar[-1][1] > MIN_CLIPPED_LENGTH:
            pos = read.reference_start if strand == -1 else read.reference_end
            read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                        read.reference_name, strand, read_length, total_segment_length, cigar[0][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None))
            read_segments[-1].is_clipped = True
            
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
            elif op_len > sv_size:
                ref_end = ref_start + ref_aligned
                read_end = read_start + read_aligned
                if read.is_reverse:
                    del_start, del_end = read_length - read_end, read_length - read_start
                else:
                    del_start, del_end = read_start , read_end
                read_segments.append(ReadSegment(align_start, del_start, del_end, ref_start, ref_end, ref_start, ref_end, read.query_name,
                                                 read.reference_name, strand, read_length,total_segment_length,read_aligned,
                                                 haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, is_primary))
                read_start = read_end+1
                ref_start = ref_end+op_len+1
                read_aligned = 0
                ref_aligned = 0

        if op == CIGAR_INS:
            if op_len < sv_size:
                read_aligned += op_len
            else:
                ins_start = read_start + read_aligned
                read_aligned += op_len
                ins_pos= ref_start + ref_aligned
                ins_end = read_start + read_aligned
                read_segments.append(ReadSegment(align_start,ins_start, ins_end, ins_pos, ins_pos, ins_pos, ins_pos, read.query_name,
                                                 read.reference_name, strand, read_length,total_segment_length,op_len, haplotype,
                                                 read.mapping_quality, genome_id, mm_rate, True, error_rate, None))
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
                                         read.mapping_quality, genome_id, mm_rate, False, error_rate, is_primary))
    merge_short_seg(read_segments)
    
        
    if is_primary and not read.has_tag('SA'):
        split = [seg for seg in read_segments if not seg.is_insertion]
        if len(split) == 1:
            read_segments = [seg for seg in read_segments if seg.is_insertion]
            read_inf = np.array([ref_ind,  read.reference_start, read.reference_end,read_length,total_segment_length, haplotype,mm_rate, error_rate, read.mapping_quality], dtype = int)
            if cigar[0][0] in CIGAR_CLIP and cigar[0][1] > MIN_CLIPPED_LENGTH:
                pos = read.reference_start if strand == 1 else read.reference_end
                read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                            read.reference_name, strand, read_length, total_segment_length, cigar[0][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None))
                read_segments[-1].is_clipped = True
            if cigar[-1][0] in CIGAR_CLIP and cigar[-1][1] > MIN_CLIPPED_LENGTH:
                pos = read.reference_start if strand == -1 else read.reference_end
                read_segments.append(ReadSegment(0, 0, 0, 0, pos, 0, pos, read.query_name,
                                            read.reference_name, strand, read_length, total_segment_length, cigar[0][1], haplotype, read.mapping_quality, genome_id, mm_rate, False, error_rate, None))
                read_segments[-1].is_clipped = True
            
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
            #if strand == -1:
            #    init_seg.read_start = seg.read_start
            #else:
            #    init_seg.read_end = seg.read_end
            #init_seg.segment_length = init_seg.read_end - init_seg.read_start
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
            read.append(ReadSegment(0, 0, s1.read_start, pos, pos, pos, pos, s1.read_id,
                                    s1.ref_id, 1, s1.read_length, s1.align_len, s1.segment_length, s1.haplotype, s1.mapq, s1.genome_id, s1.mismatch_rate, False, s1.error_rate, None))
            read[-1].is_clipped = True
        end_clip_length = s2.read_length - s2.read_end
        if end_clip_length > MIN_CLIPPED_LENGTH:
            pos = s2.ref_end if s2.strand == 1 else s2.ref_start
            read.append(ReadSegment(s2.read_end, s2.read_end, s2.read_length, pos, pos, pos, pos, s2.read_id,
                                    s2.ref_id, -1, s2.read_length, s2.align_len, s2.segment_length, s2.haplotype, s2.mapq, s2.genome_id, s2.mismatch_rate, False, s2.error_rate, None))
            read[-1].is_clipped = True
        read.sort(key=lambda s: s.read_start)
        
def get_cov(bam_file, genome_id, ref_id, poslist, min_mapq):
    region_start, region_end = max(0, poslist[0]-3), poslist[-1]+3
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    cov_list = defaultdict(list)
    BUFF = 50
    for pos in poslist:
        cov_list[(ref_id, pos)] = [0,0,0,0,0]
    for aln in aln_file.fetch(ref_id, region_start, region_end,  multiple_iterators=True):
        if not aln.is_unmapped:
            if aln.is_secondary:
                strt = bisect.bisect_left(poslist, aln.reference_start)
                end = bisect.bisect_right(poslist, aln.reference_end)
                if not strt == end:
                    for pos in poslist[strt:end]:
                        if abs(aln.reference_start - pos) <= BUFF:
                            cov_list[(ref_id, pos)][4]+=1
                            
                        elif abs(aln.reference_end - pos) <= BUFF:
                            cov_list[(ref_id, pos)][3]+=1
                            
            elif aln.mapping_quality > min_mapq:
                strt = bisect.bisect_right(poslist, aln.reference_start)
                end = bisect.bisect_left(poslist, aln.reference_end)
                if not strt == end:
                    if aln.has_tag('HP'):
                        hp = aln.get_tag('HP')
                    else:
                        hp = 0
                    for pos in poslist[strt:end]:
                        cov_list[(ref_id, pos)][hp]+=1
    return cov_list


def get_coverage_parallel(bam_files, genome_ids, thread_pool, min_mapq, double_breaks):
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
        tasks = [(bam_files[genome_id], genome_id, ref_id, pos, min_mapq) for ref_id, poslist in db_list.items() for pos in poslist]
        parsing_results = None
        parsing_results = thread_pool.starmap(get_cov, tasks)
        for item in parsing_results:
            for key, value in item.items():
                covlist[key] = value
        for db in double_breaks:
            db.bp_1.spanning_reads[genome_id] = covlist[(db.bp_1.ref_id, db.bp_1.position)]
            db.bp_2.spanning_reads[genome_id] = covlist[(db.bp_2.ref_id, db.bp_2.position)]
    
        
def get_all_reads(bam_file, region, genome_id,sv_size,use_supplementary_tag):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """

    alignments = []
    read_info_final = []
    ncol= 10000
    read_info = np.zeros((ncol,9), dtype = int)
    ref_ind, ref_id, region_start, region_end = region
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    t=0
    for aln in aln_file.fetch(ref_id, region_start, region_end,  multiple_iterators=True):
        if not aln.is_secondary and not aln.is_unmapped:
            new_segment, read_inf = get_segment(aln, genome_id, sv_size,use_supplementary_tag, ref_ind)
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


def get_all_reads_parallel(bam_file, thread_pool, ref_lengths, genome_id,
                           coverage_histograms, mismatch_histograms, n90ls, bg_mmls,read_qual,read_qual_len, args):

    CHUNK_SIZE = 10000000
    sv_size = args.sv_size
    use_supplementary_tag = args.use_supplementary_tag
    
    all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
    fetch_list = []
    for j, ctg in enumerate(all_reference_ids):
        ctg_len = ref_lengths[ctg]
        for i in range(0, max(ctg_len // CHUNK_SIZE, 1)):
            reg_start = i * CHUNK_SIZE
            reg_end = (i + 1) * CHUNK_SIZE
            if ctg_len - reg_end < CHUNK_SIZE:
                reg_end = ctg_len
            fetch_list.append((j, ctg, reg_start, reg_end))
    tasks = [(bam_file, region, genome_id,sv_size,use_supplementary_tag) for region in fetch_list]
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
    bg_mm = background_mm_rat(parsing_results)
    update_mm_hist(parsing_results, mismatch_histograms, ref_lengths)
    n90 = get_read_statistics(parsing_results, segments_by_read)
    n90 = min(n90, args.min_aligned_length)
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
                read_qual['FAIL'] += 1
                read_qual_len['FAIL'] += alignments[1][j][4]

def update_coverage_hist(coverage_histograms,genome_ids, ref_lengths, segments_by_read, control_genomes, target_genomes, loh_out):
    
    for read in segments_by_read:
        for seg in read:
            if seg.is_pass == 'PASS' and not seg.is_insertion and not seg.is_clipped:
                hist_start = seg.ref_start_ori // COV_WINDOW
                hist_end = min([seg.ref_end_ori, ref_lengths[seg.ref_id]])// COV_WINDOW
                for i in range(hist_start + 1, hist_end):
                    coverage_histograms[(seg.genome_id, seg.haplotype, seg.ref_id)][i] += 1

    for genome_id in genome_ids:
        by_hp = {}
        for hp in range(0, NUM_HAPLOTYPES):
            by_hp[hp] = []
            for chr_id, _chr_len in ref_lengths.items():
                by_hp[hp].extend(coverage_histograms[(genome_id, hp, chr_id)])

        hp1_cov, hp2_cov, hp0_cov = np.median(by_hp[1]), np.median(by_hp[2]), np.median(by_hp[0])
        logger.info(f"\tMedian coverage by PASS reads for {genome_id} (H1 / H2 / H0): {hp1_cov} / {hp2_cov} / {hp0_cov}")
        
        if loh_out:
            extract_LOH(coverage_histograms, ref_lengths, control_genomes, target_genomes, loh_out)
        

def add_read_qual(segments_by_read, ref_lengths, bg_mm, mismatch_histograms,read_qual,read_qual_len, args):
    min_mapq = args.min_mapping_quality
    min_aligned_length = args.min_aligned_length
    write_segdups_out = args.write_segdups_out
    
    mm_hist_high = background_mm_hist(segments_by_read, mismatch_histograms, bg_mm, ref_lengths)

    if write_segdups_out:
        extract_segdups(mm_hist_high, write_segdups_out) ###
    for alignments in segments_by_read:
        if not alignments:
            continue
        label_reads(alignments, min_mapq, bg_mm, mm_hist_high, min_aligned_length)

    write_readqual(segments_by_read, args.outpath_readqual,read_qual,read_qual_len)


def label_reads(read, min_mapq, bg_mm, mm_hist_high, min_aligned_length):
    MIN_ALIGNED_RATE = 0.5

    for seg in read:
        if seg.mapq < min_mapq:
            seg.is_pass += '_LOW_MAPQ'
        if high_mm_check(mm_hist_high, bg_mm, seg):
            seg.is_pass += '_HIGH_MM_rate'

    seg_ins = [1 for seg in read if not seg.is_insertion and not seg.is_clipped]
    aligned_len = sum([seg.read_end - seg.read_start for seg in read if not seg.is_clipped]) if seg_ins else read[0].align_len
    aligned_ratio = aligned_len/read[0].read_length

    if aligned_ratio < MIN_ALIGNED_RATE:
        for seg in read:
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



def background_mm_rat(parsing_results):
    COV_WINDOW_BG_MM = 2000
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
    bg_mm = np.quantile(mm_list, 0.95)
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


def extract_LOH(coverage_histograms, ref_lengths, control_genomes, target_genomes, write_loh_out):
    MIN_COV = 3
    MIN_DIFF = 2000
    MIN_LOH = 10000
    LOH_dict = defaultdict(list)
    control_genome = list(control_genomes)[0]
    for ref_id in ref_lengths.keys():
        for target_genome in list(target_genomes):
            for i, (hp1, hp2) in enumerate(zip(coverage_histograms[(target_genome, 1, ref_id)], coverage_histograms[(target_genome, 2, ref_id)])):
                if not coverage_histograms[(control_genome, 1, ref_id)][i] > MIN_COV and coverage_histograms[(control_genome, 2, ref_id)][i] > MIN_COV:
                    continue
                if (hp1 <= MIN_COV and hp2 > MIN_COV):
                    LOH_dict[(target_genome, 1, ref_id)].append(i)
                elif (hp1 > MIN_COV and hp2 <= MIN_COV) :
                    LOH_dict[(target_genome, 2, ref_id)].append(i)
    LOH_region = defaultdict(list)
    for key, loh_reg in LOH_dict.items():
        for ind in loh_reg:
            if not LOH_region[key]:
                LOH_region[key].append([ind * COV_WINDOW])
                LOH_region[key].append([(ind + 1) * COV_WINDOW])
            elif ind * COV_WINDOW - LOH_region[key][1][-1] <= MIN_DIFF:
                LOH_region[key][1][-1]=(ind + 1) * COV_WINDOW
            else:
                LOH_region[key][0].append(ind * COV_WINDOW)
                LOH_region[key][1].append((ind + 1) * COV_WINDOW)
    
    write_loh_out.write("#chr_id\tstart\tend\tgenome_ids\thaplotype\t\n")  
    for key, loh_reg in LOH_region.items():
        for (st,end) in zip(loh_reg[0], loh_reg[1]):
            if end - st >= MIN_LOH:
                write_loh_out.write('\t'.join([key[2], str(st), str(end), key[0], str(key[1])]) + '\n')
            
 
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
            if aln.is_clipped or aln.is_insertion:
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
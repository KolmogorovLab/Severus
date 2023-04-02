import pysam
import numpy as np
from collections import  defaultdict
import time
import bisect

class ReadSegment(object):
    __slots__ = ("read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                 "strand", "read_length",'segment_length', "haplotype", "mapq", "genome_id",
                 'mismatch_rate', 'is_pass', "is_insertion", "is_clipped", 'is_end', 'bg_mm_rate')
    def __init__(self, read_start, read_end, ref_start, ref_end, read_id, ref_id,
                 strand, read_length,segment_length, haplotype, mapq, genome_id, mismatch_rate, is_insertion):
        self.read_start = read_start
        self.read_end = read_end
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_id = read_id
        self.ref_id = ref_id
        self.strand = strand
        self.read_length = read_length
        self.segment_length = segment_length
        self.haplotype = haplotype
        self.mapq = mapq
        self.genome_id = genome_id
        self.mismatch_rate = mismatch_rate
        self.is_pass = 'PASS'
        self.is_insertion = is_insertion
        self.is_clipped = False
        self.is_end = False
        self.bg_mm_rate = 1
    def __str__(self):
        return "".join(["read_start=", str(self.read_start), " read_end=", str(self.read_end), " ref_start=", str(self.ref_start),
                         " ref_end=", str(self.ref_end), " read_id=", str(self.read_id), " ref_id=", str(self.ref_id), " strand=", str(self.strand),
                         " read_length=", str(self.read_length), " haplotype=", str(self.haplotype),
                         " mapq=", str(self.mapq), "mismatch_rate=", str(self.mismatch_rate), " read_qual=", str(self.is_pass), "genome_id=", str(self.genome_id)])
    
def get_segment(read, genome_id,sv_size):
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
    
    ref_start = read.reference_start
    read_segments =[]
    cigar = read.cigartuples
    
    
    num_of_mismatch = 0
    nm = read.get_tag('NM')
    indel = sum([b for a, b in cigar if a in [CIGAR_INS, CIGAR_DEL]])
    num_of_mismatch = nm - indel 
    total_segment_length = sum([b for a, b in cigar if a not in [CIGAR_CLIP, CIGAR_DEL]])
    mm_rate = num_of_mismatch / total_segment_length
    #mm_rate = read.get_tag('de')
    read_length = np.sum([k[1] for k in cigar if k[0] != CIGAR_DEL])
    strand = '-' if read.is_reverse else '+'
    
    if read.has_tag('HP'):
        haplotype = read.get_tag('HP')
    else:
        haplotype = 0
        
    for token in cigar:
        op = token[0]
        op_len = token[1]
        if op in CIGAR_CLIP:
            if first_clip:
                read_start = op_len
        first_clip = False
        if op in CIGAR_MATCH:
            read_aligned += op_len
            ref_aligned += op_len
            #if op == CIGAR_MM:
                #num_of_mismatch += op_len
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
                #mm_rate = num_of_mismatch / read_aligned
                read_segments.append(ReadSegment(del_start, del_end, ref_start, ref_end, read.query_name,
                                        read.reference_name, strand, read_length,read_aligned, haplotype, read.mapping_quality, genome_id, mm_rate, False))
                read_start = read_end+1
                ref_start = ref_end+op_len+1
                read_aligned = 0
                ref_aligned = 0
                #num_of_mismatch = 0
        if op == CIGAR_INS:
            if op_len < sv_size:
                read_aligned += op_len
            else:
                ins_start = read_start +read_aligned
                read_aligned += op_len
                ins_pos= ref_start + ref_aligned
                ins_end = read_start +read_aligned
                #mm_rate = 0
                read_segments.append(ReadSegment(ins_start, ins_end, ins_pos, ins_pos, read.query_name,
                                        read.reference_name, strand, read_length,op_len, haplotype, read.mapping_quality, genome_id,mm_rate, True))
    if ref_aligned !=0:
        ref_end = ref_start + ref_aligned
        read_end = read_start + read_aligned
        #mm_rate = num_of_mismatch / read_aligned
        if read.is_reverse:
            read_start, read_end = read_length - read_end, read_length - read_start
        read_segments.append(ReadSegment(read_start, read_end, ref_start, ref_end, read.query_name,
                       read.reference_name, strand, read_length,read_aligned, haplotype, read.mapping_quality, genome_id,mm_rate, False))
    return read_segments 


def add_clipped_end(segments_by_read):
    MAX_CLIPPED_LENGTH = 200
    
    for read in segments_by_read.values():
        read2 = [seg for seg  in read if not seg.is_insertion]
        read2.sort(key=lambda s: s.read_start)
        s1 = read2[0]
        s2 = read2[-1]
        s1.is_end = True
        s2.is_end = True
        if s1.read_start > MAX_CLIPPED_LENGTH:
            pos = s1.ref_start if s1.strand == '+' else s1.ref_end
            read.append(ReadSegment(0, s1.read_start, pos, pos, s1.read_id,
                                    s1.ref_id, s1.strand, s1.read_length,0, s1.haplotype, s1.mapq, s1.genome_id,0, False))
            read[-1].is_clipped = True
        end_clip_length = s2.read_length - s2.read_end
        if end_clip_length > MAX_CLIPPED_LENGTH:
            pos = s2.ref_end if s2.strand == '+' else s2.ref_start
            read.append(ReadSegment(s2.read_end, s2.read_length, pos, pos, s2.read_id,
                                    s2.ref_id, s2.strand, s2.read_length, 0, s2.haplotype, s2.mapq, s2.genome_id,0, False))
            read[-1].is_clipped = True
        read.sort(key=lambda s: s.read_start)
    return segments_by_read
            

def get_all_reads(bam_file, region, genome_id,sv_size):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """
    alignments = []
    ref_id, region_start, region_end = region
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(ref_id, region_start, region_end,  multiple_iterators=True):
        if not aln.is_secondary and not aln.is_unmapped:
            new_segment= get_segment(aln, genome_id, sv_size)
            for new_seg in new_segment:
                alignments.append(new_seg)
    return alignments


def get_all_reads_parallel(bam_file, thread_pool, ref_lengths, genome_id,
                           min_mapq, sv_size):
    CHUNK_SIZE = 10000000
    all_reference_ids = [r for r in pysam.AlignmentFile(bam_file, "rb").references]
    fetch_list = []
    for ctg in all_reference_ids:
        ctg_len = ref_lengths[ctg]
        for i in range(0, max(ctg_len // CHUNK_SIZE, 1)):
            reg_start = i * CHUNK_SIZE
            reg_end = (i + 1) * CHUNK_SIZE
            if ctg_len - reg_end < CHUNK_SIZE:
                reg_end = ctg_len
            fetch_list.append((ctg, reg_start, reg_end))
            
    tasks = [(bam_file, region, genome_id,sv_size) for region in fetch_list]
    parsing_results = None
    parsing_results = thread_pool.starmap(get_all_reads, tasks)
    segments_by_read = defaultdict(list)
    for alignments in parsing_results:
        for aln in alignments:
            segments_by_read[aln.read_id].append(aln)
    segments_by_read = add_clipped_end(segments_by_read)
    return segments_by_read
        
           
COV_WINDOW = 500
def update_coverage_hist(genome_ids, ref_lengths, segments_by_read):
    NUM_HAPLOTYPES = 3
    coverage_histograms = {}
    
    for genome_id in genome_ids:
        for chr_id, chr_len in ref_lengths.items():
            for hp in range(0,NUM_HAPLOTYPES):
                coverage_histograms[(genome_id, hp, chr_id)] = [0 for _ in range(chr_len // COV_WINDOW + 1)]
                
    for read in segments_by_read.values():
        for seg in read:
            if seg.is_pass == 'PASS' and not seg.is_insertion and not seg.is_clipped:
                hist_start = seg.ref_start // COV_WINDOW
                hist_end = seg.ref_end // COV_WINDOW
                for i in range(hist_start + 1, hist_end): ## Check with (hist_start, hist_end + 1)
                    coverage_histograms[(seg.genome_id, seg.haplotype, seg.ref_id)][i] += 1
                    
    return coverage_histograms


def filter_reads(segments_by_read):
    segments_by_read_filtered = []
    for read in segments_by_read.values():
        new_read = []
        for seg in read:
            if seg.is_pass == 'PASS':
                new_read.append(seg)
        segments_by_read_filtered.append(new_read)
    return segments_by_read_filtered
        

def add_read_qual(segments_by_read, ref_lengths, thread_pool, min_mapq, max_error_rate):
    mm_hist_low, mm_hist_high = background_mm_hist(segments_by_read, ref_lengths)
    high_mm_region = extract_segdups(mm_hist_low, mm_hist_high)
    for read, alignments in segments_by_read.items():
        segments_by_read[read] = label_reads(alignments, min_mapq, mm_hist_low, high_mm_region, max_error_rate)
    
    
def label_reads(read, min_mapq, mm_hist_low, high_mm_region, max_error_rate):
    
    MIN_ALIGNED_RATE = 0.5
    MIN_ALIGNED_LENGTH = 5000
    MAX_SEGMENTED_READ = 10
    MAX_CHR_SPAN = 2
    MIN_SEGMENT_LEN = 100
    #MAX_MM_RATE_THR = 0.005
    #MAX_ERROR_RATE = 0.01
    fail_type = ''
    chr_list = []
    dedup_segments = []
    dedup_segments_ins = []
    n_seg = 0#
    
    for seg in read:
        if not dedup_segments or dedup_segments[-1].read_start != seg.read_start:            
            if seg.is_insertion or seg.is_clipped:
                dedup_segments_ins.append(seg)
            else:
                dedup_segments.append(seg)
                chr_list.append(seg.ref_id)
                n_seg += 1  
    dedup_segments +=dedup_segments_ins        
    aligned_len = sum([seg.segment_length for seg in dedup_segments if not seg.is_insertion and not seg.is_clipped])
    aligned_ratio = aligned_len/read[0].read_length
    
    if n_seg > MAX_SEGMENTED_READ or len(set(chr_list)) > MAX_CHR_SPAN:
    #if n_seg > MAX_SEGMENTED_READ:
        fail_type = 'SEGMENTED'
    elif aligned_len < MIN_ALIGNED_LENGTH or aligned_ratio < MIN_ALIGNED_RATE:
        fail_type = 'LOW_ALIGNED_LEN'#
        
    if fail_type:
        for seg in dedup_segments:
            seg.is_pass = fail_type
        return dedup_segments#
    
    for seg in dedup_segments:
        if not seg.is_insertion and not seg.is_clipped and seg.segment_length < MIN_SEGMENT_LEN:
            seg.is_pass = 'SEG_LENGTH'
            continue
        if seg.mapq < min_mapq:
            seg.is_pass = 'LOW_MAPQ'
            continue
        high_mm_check(high_mm_region, mm_hist_low, seg)
        if seg.mismatch_rate > seg.bg_mm_rate + max_error_rate:
            seg.is_pass = 'HIGH_MM_rate'
            continue#
            
    return dedup_segments


def get_background_mm_rate(seg, mismatch_histograms):
    hist_start = seg.ref_start // COV_WINDOW
    hist_end = seg.ref_end // COV_WINDOW
    mm_list = mismatch_histograms[seg.ref_id][hist_start : hist_end + 1]
    if mm_list:
        seg.bg_mm_rate = float(np.median(mm_list))
    else:
        seg.bg_mm_rate = 1
        

def background_mm_hist(segments_by_read, ref_lengths):
    COV_WINDOW  = 500
    MED_THR = 3
    mismatch_histograms = defaultdict(list)
    mm_hist_low = {}
    mm_hist_high = {}
    
    for chr_id, chr_len in ref_lengths.items():
        mismatch_histograms[chr_id] = [[] for _ in range(chr_len // COV_WINDOW + 2)]
        mm_hist_low[chr_id] = [1 for _ in range(chr_len // COV_WINDOW + 2)]
        mm_hist_high[chr_id] = [1 for _ in range(chr_len // COV_WINDOW + 2)]
        
    for read in segments_by_read.values():
        for seg in read:
            if seg.is_clipped or seg.is_insertion:
                continue
            hist_start = seg.ref_start // COV_WINDOW
            hist_end = seg.ref_end // COV_WINDOW
            for i in range(hist_start, hist_end + 1):
                mismatch_histograms[seg.ref_id][i].append(seg.mismatch_rate)
                
    for chr_id, mm_rate in mismatch_histograms.items():
        for i, mm_list in enumerate(mm_rate):
            if not mm_list or len(mm_list)<5:
                mm_hist_low[chr_id][i] = 1
                mm_hist_high[chr_id][i] = 1
                continue
            mm_list.sort()
            mm_hist_low[chr_id][i] = mm_list[MED_THR]
            mm_hist_high[chr_id][i] = mm_list[-MED_THR]
    return mm_hist_low, mm_hist_high


def extract_segdups(mm_hist_low, mm_hist_high):
    THR = 0.005
    COV_WINDOW = 500
    high_mm_region = defaultdict(list)
    for key, mm_low_chrom in mm_hist_low.items():
        mm_high_chrom = mm_hist_high[key]
        for ind, (a, b)  in enumerate(zip(mm_low_chrom, mm_high_chrom)):
            if b - a >= THR:
                if not high_mm_region[key]:
                    high_mm_region[key].append([ind * COV_WINDOW])
                    high_mm_region[key].append([(ind + 1) * COV_WINDOW])
                elif ind * COV_WINDOW == high_mm_region[key][-1]:
                    high_mm_region[key][1][-1]=(ind + 1) * COV_WINDOW
                else:
                    high_mm_region[key][0].append(ind * COV_WINDOW)
                    high_mm_region[key][1].append((ind + 1) * COV_WINDOW)
    return high_mm_region
                    

def high_mm_check(high_mm_region, mm_hist_low, seg):
    high_mm_reg = high_mm_region[seg.ref_id]
    if high_mm_reg:
        strt = bisect.bisect_right(high_mm_reg[0], seg.ref_start)
        strt_2  = bisect.bisect_right(high_mm_reg[1], seg.ref_start)
        end = bisect.bisect_left(high_mm_reg[0],seg.ref_end)
        end_1 = bisect.bisect_left(high_mm_reg[1],seg.ref_end)
        if not strt == strt_2 or not end == end_1:
            get_background_mm_rate(seg, mm_hist_low)
        
  

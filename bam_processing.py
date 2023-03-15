import pysam
from multiprocessing import Pool
import numpy as np
from collections import  defaultdict

class ReadSegment(object):
    __slots__ = ("read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                 "strand", "read_length",'segment_length', "haplotype", "mapq", "genome_id",
                 'mismatch_rate', "is_insertion", "is_clipped", 'is_end')
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
        self.is_insertion = is_insertion
        self.is_clipped = False
        self.is_end = False
    def __str__(self):
        return "".join(["read_start=", str(self.read_start), " read_end=", str(self.read_end), " ref_start=", str(self.ref_start),
                         " ref_end=", str(self.ref_end), " read_id=", str(self.read_id), " ref_id=", str(self.ref_id), " strand=", str(self.strand),
                         " read_length=", str(self.read_length), " haplotype=", str(self.haplotype),
                         " mapq=", str(self.mapq), " genome_id=", str(self.genome_id)])
    
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
    
    #num_of_mismatch = 0
    nm = read.get_tag('NM')
    indel = sum([b for a, b in cigar if a in [CIGAR_INS, CIGAR_DEL] and b < sv_size])
    num_of_mismatch = nm - indel 
    total_segment_length = sum([b for a, b in cigar if a not in [CIGAR_CLIP, CIGAR_DEL]])
    mm_rate = num_of_mismatch / total_segment_length
    
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
                num_of_mismatch = 0
        if op == CIGAR_INS:
            if op_len < sv_size:
                read_aligned += op_len
            else:
                ins_start = read_start +read_aligned
                read_aligned += op_len
                ins_pos= ref_start + ref_aligned
                ins_end = read_start +read_aligned
                mm_rate = 0
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
                           min_mapq, sv_size, add_clipped_ends):
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
    #thread_pool = Pool(num_threads)
    parsing_results = thread_pool.starmap(get_all_reads, tasks)
    segments_by_read = defaultdict(list)
    for alignments in parsing_results:
        for aln in alignments:
            segments_by_read[aln.read_id].append(aln)
    if add_clipped_ends:
        segments_by_read = add_clipped_end(segments_by_read)
    return segments_by_read        
       
            
COV_WINDOW = 500
def update_coverage_hist(genome_ids, ref_lengths, segments_by_read, min_mapq, max_read_error):
    NUM_HAPLOTYPES = 3
    segments_by_read_filtered = filter_all_reads(segments_by_read, min_mapq, max_read_error)
    allsegments = get_allsegments(segments_by_read_filtered)
    coverage_histograms = {}
    for genome_id in genome_ids:
        for chr_id, chr_len in ref_lengths.items():
            for hp in range(0,NUM_HAPLOTYPES):
                coverage_histograms[(genome_id, hp, chr_id)] = [0 for _ in range(chr_len // COV_WINDOW + 1)]
    
    for read in allsegments:
        hist_start = read.ref_start // COV_WINDOW
        hist_end = read.ref_end // COV_WINDOW
        for i in range(hist_start + 1, hist_end): ## Check with (hist_start, hist_end + 1)
            coverage_histograms[(read.genome_id, read.haplotype, read.ref_id)][i] += 1
    return coverage_histograms


def filter_all_reads(segments_by_read,min_mapq, max_read_error):
    MIN_ALIGNED_LENGTH = 5000
    MIN_ALIGNED_RATE = 0.5
    MAX_SEGMENTS = 10
    MIN_SEGMENT_LENGTH = 100
    
    segments_by_read_filtered=[]
    for read_id, segments in segments_by_read.items():
        dedup_segments = []
        segments.sort(key=lambda s: s.read_start)
        
        for seg in segments:
            if not dedup_segments or dedup_segments[-1].read_start != seg.read_start:            
                if (seg.is_insertion or seg.is_clipped) and seg.mapq > min_mapq:
                    dedup_segments.append(seg)
                elif seg.mapq>min_mapq and seg.segment_length > MIN_SEGMENT_LENGTH and seg.mismatch_rate < max_read_error:
                    dedup_segments.append(seg)
                    
        aligned_len = sum([seg.segment_length for seg in dedup_segments if not seg.is_insertion and not seg.is_clipped])
        aligned_ratio = aligned_len/segments[0].read_length
        if aligned_len < MIN_ALIGNED_LENGTH or aligned_ratio < MIN_ALIGNED_RATE or len(segments) > MAX_SEGMENTS:
            continue
        segments_by_read_filtered.append(dedup_segments)
        
    return segments_by_read_filtered


def get_allsegments(segments_by_read_filtered):
    allsegments = []
    for read in segments_by_read_filtered:
        for seg in read:
            if not seg.is_insertion:
                allsegments.append(seg)
    return allsegments
    





    

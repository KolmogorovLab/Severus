import pysam
from multiprocessing import Pool
import numpy as np
from collections import  defaultdict



class ReadSegment(object):
    __slots__ = ("read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                 "strand", "read_length",'segment_length', "haplotype", "mapq", "genome_id",'mismatch_rate', "is_insertion")
    def __init__(self, read_start, read_end, ref_start, ref_end, read_id, ref_id,
                 strand, read_length,segment_length, haplotype, mapq, genome_id,mismatch_rate, is_insertion):
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
    def __str__(self):
        return "".join(["read_start=", str(self.read_start), " read_end=", str(self.read_end), " ref_start=", str(self.ref_start),
                         " ref_end=", str(self.ref_end), " read_id=", str(self.read_id), " ref_id=", str(self.ref_id), " strand=", str(self.strand),
                         " read_length=", str(self.read_length), " haplotype=", str(self.haplotype),
                         " mapq=", str(self.mapq), " genome_id=", str(self.genome_id)])
    
def get_segment(read, genome_id,sv_size):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    ref_start = read.reference_start
    read_segments =[]
    cigar = read.cigartuples
    read_length = np.sum([k[1] for k in cigar if k[0]!=2])
    num_of_mismatch = 0
    strand = '-' if read.is_reverse else '+'
    if read.has_tag('HP'):
        haplotype = read.get_tag('HP')
    else:
        haplotype = 0
    for token in cigar:
        op = token[0]
        op_len = token[1]
        if op in [4, 5]:
            if first_clip:
                read_start = op_len
        first_clip = False
        if op in [0,7,8]:
            read_aligned += op_len
            ref_aligned += op_len
        if op == 8:
            num_of_mismatch +=1
        if op == 2:
            if op_len < sv_size:
                ref_aligned += op_len
            elif op_len > sv_size:
                ref_end = ref_start + ref_aligned
                read_end = read_start + read_aligned
                if read.is_reverse:
                    del_start, del_end = read_length - read_end, read_length - read_start
                else:
                    del_start, del_end = read_start , read_end
                mm_rate = num_of_mismatch/read_aligned
                read_segments.append(ReadSegment(del_start, del_end, ref_start, ref_end, read.query_name,
                                        read.reference_name, strand, read_length,read_aligned, haplotype, read.mapping_quality, genome_id, mm_rate,False))
                read_start = read_end+1
                ref_start = ref_end+op_len+1
                read_aligned = 0
                ref_aligned = 0
                num_of_mismatch = 0
        if op == 1:
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
        mm_rate = num_of_mismatch/read_aligned
        if read.is_reverse:
            read_start, read_end = read_length - read_end, read_length - read_start
        read_segments.append(ReadSegment(read_start, read_end, ref_start, ref_end, read.query_name,
                       read.reference_name, strand, read_length,read_aligned, haplotype, read.mapping_quality, genome_id,mm_rate, False))
    return read_segments 


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


def get_all_reads_parallel(bam_file, thread_pool, ref_lengths,
                           min_mapq, genome_id,sv_size):
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
    return segments_by_read

def extract_lowmapq_regions(segments_by_read , min_mapq):
    MAX_LOWMAPQ_READS = 10
    WIN_SIZE = 2000
        
    lowmapq_reads = defaultdict(int)
    for segments in segments_by_read.values():
        segments.sort(key=lambda s: s.read_start)
        for seg in segments:
            if seg.mapq < min_mapq and not seg.is_insertion:
                for i in range(seg.ref_start//WIN_SIZE,seg.ref_end//WIN_SIZE): ###Insertion
                    lowmapq_reads[(seg.ref_id,i)]+=1
                    
    lowmapq_reads = dict(sorted(lowmapq_reads.items()))
    lowmapq_reg = defaultdict(list)
    for key, value in lowmapq_reads.items():
        if value > MAX_LOWMAPQ_READS:
            if not lowmapq_reg[key[0]]:
                lowmapq_reg[key[0]].append([key[1]*WIN_SIZE])
                lowmapq_reg[key[0]].append([(key[1]+1)*WIN_SIZE])
            elif (key[1]*WIN_SIZE)==lowmapq_reg[key[0]][1][-1]:
                lowmapq_reg[key[0]][1][-1]=(key[1]+1)*WIN_SIZE
            else:
                lowmapq_reg[key[0]][0].append(key[1]*WIN_SIZE)
                lowmapq_reg[key[0]][1].append((key[1]+1)*WIN_SIZE)
    return lowmapq_reg

        
def filter_allreads(segments_by_read,min_mapq, max_read_error):
    MIN_ALIGNED_LENGTH = 5000
    MIN_ALIGNED_RATE = 0.5
    MAX_SEGMENTS = 10
    MIN_SEGMENT_LENGTH = 100
    
    segments_by_read_filtered=[]
    for read_id, segments in segments_by_read.items():
        dedup_segments = []
        segments.sort(key=lambda s: s.read_start)
        aligned_len = sum([seg.segment_length for seg in segments if not seg.is_insertion])
        aligned_ratio = aligned_len/segments[0].read_length
        if aligned_len < MIN_ALIGNED_LENGTH or aligned_ratio < MIN_ALIGNED_RATE or len(segments) > MAX_SEGMENTS:
            continue
        
        for seg in segments:
            if not dedup_segments or dedup_segments[-1].read_start != seg.read_start:            
                if seg.is_insertion and seg.mapq>min_mapq:
                    dedup_segments.append(seg)
                elif seg.mapq>min_mapq and seg.segment_length > MIN_SEGMENT_LENGTH and seg.mismatch_rate< max_read_error:
                    dedup_segments.append(seg)
        segments_by_read_filtered.append(dedup_segments)
    return segments_by_read_filtered


def all_reads_position(allsegments):
    allreads_pos = defaultdict(list)
    for seg in allsegments:
        if not allreads_pos[seg.ref_id]:
            allreads_pos[seg.ref_id].append([seg.ref_start])
            allreads_pos[seg.ref_id].append([seg.ref_end])
            allreads_pos[seg.ref_id].append([seg.read_id])
            allreads_pos[seg.ref_id].append([(seg.haplotype,seg.genome_id)])
            allreads_pos[seg.ref_id].append([(seg.read_id)])
        else:    
            allreads_pos[seg.ref_id][0].append(seg.ref_start)
            allreads_pos[seg.ref_id][1].append(seg.ref_end)
            allreads_pos[seg.ref_id][2].append(seg.read_id)
            allreads_pos[seg.ref_id][3].append((seg.haplotype,seg.genome_id))
            allreads_pos[seg.ref_id][4].append((seg.read_id))
    return allreads_pos


def get_insertionreads(segments_by_read_filtered):
    ins_list_all = []
    for read in segments_by_read_filtered:
        for seg in read:
            if seg.is_insertion:
                ins_list_all.append(seg)
    return ins_list_all
    

def get_allsegments(segments_by_read_filtered):
    allsegments = []
    for read in segments_by_read_filtered:
        for seg in read:
            if not seg.is_insertion:
                allsegments.append(seg)
    return allsegments


def get_splitreads(segments_by_read_filtered):
    split_reads = []
    for read in segments_by_read_filtered:
        split_reads_add = []
        if len(read)>1:
            for seg in read:
                if not seg.is_insertion:
                    split_reads_add.append(seg)
            split_reads.append(split_reads_add)
    return split_reads


def resolve_overlaps(split_reads, min_ovlp_len):
    """
    Some supplementary alignments may be overlapping (e.g. in case of inversions with flanking repeat).
    This function checks if the overlap has ok structe, trims and outputs non-overlapping alignments
    """
    def _get_ovlp(seg_1, seg_2):
        max_ovlp_len = min(seg_1.read_end - seg_1.read_start, seg_2.read_end - seg_2.read_start)
        if (seg_1.read_end - seg_2.read_start > min_ovlp_len and
            seg_1.read_end - seg_2.read_start < max_ovlp_len and
            seg_2.read_end > seg_1.read_end and
            seg_2.read_start > seg_1.read_start):
            return seg_1.read_end - seg_2.read_start
        else:
            return 0
        
    new_reads = []
    for read_segments in split_reads:
        upd_segments = []
        for i in range(len(read_segments)):
            left_ovlp = 0
            if i > 0 and read_segments[i - 1].ref_id == read_segments[i].ref_id:
                left_ovlp = _get_ovlp(read_segments[i - 1], read_segments[i])
            left_ovlp = left_ovlp
            seg = read_segments[i]
            if left_ovlp > 0:
                if seg.strand == "+":
                    seg.read_start = seg.read_start + left_ovlp
                    seg.ref_start = seg.ref_start + left_ovlp
                else:
                    seg.read_start = seg.read_start + left_ovlp
                    seg.ref_end = seg.ref_end - left_ovlp
            upd_segments.append(seg)
        new_reads.append(upd_segments)
    return new_reads    
            
            
COV_WINDOW = 500
def update_coverage_hist(genome_ids, ref_lengths, NUM_HAPLOTYPES, read_alignments):
    coverage_histograms = {}
    for genome_id in genome_ids:
        for chr_id, chr_len in ref_lengths.items():
            for hp in range(0,NUM_HAPLOTYPES):
                coverage_histograms[(genome_id, hp, chr_id)] = [0 for _ in range(chr_len // COV_WINDOW + 1)]
    
    for read in read_alignments:
        hist_start = read.ref_start // COV_WINDOW
        hist_end = read.ref_end // COV_WINDOW
        for i in range(hist_start, hist_end + 1):
            coverage_histograms[(read.genome_id, read.haplotype, read.ref_id)][i] += 1
    return coverage_histograms


def segment_coverage(histograms, genome_id, ref_id, ref_start, ref_end, haplotype):
    hist_start = ref_start // COV_WINDOW
    hist_end = ref_end // COV_WINDOW
    cov_list = histograms[(genome_id, haplotype, ref_id)][hist_start : hist_end + 1]

    if not cov_list:
        return 0
    return int(np.median(cov_list))
    
def write_alignments(allsegments, outpath):
    aln_dump_stream = open(outpath, "w")
    for read in allsegments:
        aln_dump_stream.write(str(read) + "\n")
        for seg in read:
            aln_dump_stream.write(str(seg) + "\n")
    aln_dump_stream.write("\n")



    





    

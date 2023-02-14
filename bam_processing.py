import sys
import re
import pysam
from multiprocessing import Pool
import random
import os
import subprocess
import numpy as np
from collections import namedtuple, defaultdict,Counter

#read alignment constants
MIN_ALIGNED_LENGTH = 5000
MIN_ALIGNED_RATE = 0.5
MAX_SEGMENTS = 10
MIN_SEGMENT_LENGTH = 100
min_lowmapq_reads = 10


#ReadSegment = namedtuple("ReadSegment", ["read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
#                                         "strand", "read_length", "haplotype", "mapq", "genome_id"])
class ReadSegment(object):
    __slots__ = ("read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                 "strand", "read_length", "haplotype", "mapq", "genome_id")
    def __init__(self, read_start, read_end, ref_start, ref_end, read_id, ref_id,
                 strand, read_length, haplotype, mapq, genome_id):
        self.read_start = read_start
        self.read_end = read_end
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.read_id = read_id
        self.ref_id = ref_id
        self.strand = strand
        self.read_length = read_length
        self.haplotype = haplotype
        self.mapq = mapq
        self.genome_id = genome_id
    def __str__(self):
        return "".join(["read_start=", str(self.read_start), " read_end=", str(self.read_end), " ref_start=", str(self.ref_start),
                         " ref_end=", str(self.ref_end), " read_id=", str(self.read_id), " ref_id=", str(self.ref_id), " strand=", str(self.strand),
                         " read_length=", str(self.read_length), " haplotype=", str(self.haplotype),
                         " mapq=", str(self.mapq), " genome_id=", str(self.genome_id)])



#TODO: merge with get_segment below
ReadSegmentLegacy = namedtuple("ReadSegmentLegacy", ["read_start", "read_end", "ref_start", "ref_end", "read_id", "ref_id",
                                                     "strand", "read_length", "mismatch_rate", 'mapq'])
def _get_segment_legacy(read_id, ref_id, ref_start, strand, cigar, num_mismatch, mapq):
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    for token in cigar_parser.findall(cigar):
        op = token[-1]
        op_len = int(token[:-1])

        if op == "H" or op == "S":
            if first_clip:
                read_start = op_len
            read_length += op_len
        first_clip = False
        if op == "M" or op == "=" or op == "X":
            read_aligned += op_len
            ref_aligned += op_len
            read_length += op_len
        if op == 'D':
            ref_aligned += op_len
        if op == "I":
            read_aligned += op_len
            read_length += op_len
    ref_end = ref_start + ref_aligned
    read_end = read_start + read_aligned
    length_diff = abs(ref_aligned - read_aligned)
    mismatch_rate = (num_mismatch - length_diff) / (read_aligned + 1)
    if strand == "-":
        read_start, read_end = read_length - read_end, read_length - read_start
    return ReadSegmentLegacy(read_start, read_end, ref_start, ref_end, read_id,
                             ref_id, strand, read_length, mismatch_rate, mapq)
###

def check_read_mapping_confidence(read_id, flags, ref_id, ref_start,strand,cigar, mapq, num_mismatches, sa_tag,min_aln_length, min_aligned_rate,
                                  max_read_error, max_segments):
    read_info = []
    is_secondary = int(flags) & 0x100
    segments = [_get_segment_legacy(read_id, ref_id, ref_start, strand, cigar, num_mismatches,mapq)]
    if sa_tag:
        for sa_aln in sa_tag.split(";"):
            if sa_aln:
                sa_fields = sa_aln.split(",")
                sa_ref, sa_ref_pos, sa_strand, sa_cigar, sa_mismatches, sa_mapq= \
                        sa_fields[0], int(sa_fields[1]), sa_fields[2], sa_fields[3], int(sa_fields[5]),int(sa_fields[4])
                segments.append(_get_segment_legacy(read_id, sa_ref, sa_ref_pos, sa_strand, sa_cigar, sa_mismatches, sa_mapq))
    segments.sort(key=lambda s: s.read_start)
    read_length = segments[0].read_length
    WND_LEN = 100
    read_coverage = [0 for x in range(read_length // WND_LEN)]
    weighted_mm_sum = 0
    total_segment_length = 0
    for seg in segments:
        for i in range(seg.read_start // WND_LEN, seg.read_end // WND_LEN):
            read_coverage[i] = 1
        weighted_mm_sum += seg.mismatch_rate * (seg.read_end - seg.read_start)
        total_segment_length += seg.read_end - seg.read_start
        read_info.append(seg)
    aligned_length = sum(read_coverage) * WND_LEN
    mean_mm = weighted_mm_sum / (total_segment_length + 1)
    if (is_secondary or 
            aligned_length < min_aln_length or 
            aligned_length / read_length < min_aligned_rate or
            len(segments) > max_segments or
            mean_mm > max_read_error):
        return False , read_info
    else:
        return True, read_info


cigar_parser = re.compile("[0-9]+[MIDNSHP=X]")
def get_segment(read_id, ref_id, ref_start, strand, cigar, haplotype, mapq, genome_id,sv_size):
    """
    Parses cigar and generate ReadSegment structure with alignment coordinates
    """
    first_clip = True
    read_start = 0
    read_aligned = 0
    read_length = 0
    ref_aligned = 0
    #read_end = 0
    read =[]
    ins_list=[]
    add_del = []

    for token in cigar_parser.findall(cigar):
        op = token[-1]
        op_len = int(token[:-1])

        if op == "H" or op == "S":
            if first_clip:
                read_start = op_len
            read_length += op_len
        first_clip = False

        if op == "M" or op == "=" or op == "X":
            read_aligned += op_len
            ref_aligned += op_len
            read_length += op_len
        if op == 'D':
            if op_len < sv_size:
                ref_aligned += op_len
            elif op_len > sv_size:
                ref_end = ref_start + ref_aligned
                read_end = read_start + read_aligned
                add_del.append((read_start, read_end,ref_start, ref_end))
                read_start = read_end+1
                read_aligned = 0
                ref_aligned = 0
                ref_start = ref_end+op_len
        if op == 'I':
            if op_len < sv_size:
                read_aligned += op_len
                read_length += op_len                
            else:
                read_aligned += op_len
                read_length += op_len
                ref_end1 = ref_start + ref_aligned
                ins_list.append([ref_id,ref_end1,read_id,genome_id,haplotype,op_len])
    if ref_aligned !=0:
        ref_end = ref_start + ref_aligned
        read_end = read_start + read_aligned

        if strand == "-":
            read_start, read_end = read_length - read_end, read_length - read_start
        read.append(ReadSegment(read_start, read_end, ref_start, ref_end, read_id,
                       ref_id, strand, read_length, haplotype, mapq, genome_id))
        if add_del:
            if strand == "-":
                for seg in add_del:
                    read_start, read_end = read_length - seg[1], read_length - seg[0]
                    read.append(ReadSegment(read_start, read_end, seg[2], seg[3], read_id,
                                            ref_id, strand, read_length, haplotype, mapq, genome_id))
            else:
                for seg in add_del:
                    read.append(ReadSegment(seg[0], seg[1], seg[2], seg[3], read_id,
                               ref_id, strand, read_length, haplotype, mapq, genome_id))
    return read, ins_list




def get_split_reads(bam_file, region, max_read_error, min_mapq, genome_id,sv_size):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """
    all_reads = []
    alignments = []
    lowmapq = []
    ref_id, region_start, region_end = region
    inslist = []
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(ref_id, region_start, region_end,  multiple_iterators=True):
        
        fields = aln.to_string().split()
        read_id, flags, ref_id = fields[0:3]
        cigar = fields[5]
        mapq = int(fields[4])
        ref_start = int(fields[3])
        is_supplementary = int(flags) & 0x800
        is_secondary = int(flags) & 0x100
        is_unmapped = int(flags) & 0x4
        is_primary = not (is_supplementary or is_secondary or is_unmapped)
        strand = "-" if int(flags) & 0x10 else "+"
        sa_tag = None
        hp_tag = None
        for tag in fields[11:]:
            if tag.startswith("SA"):
                sa_tag = tag[5:]
            if tag.startswith("HP"):
                hp_tag = tag[5:]
            if tag.startswith("NM"):
                num_mismatches = int(tag[5:])
        if not hp_tag:
            haplotype = 0
        else:
            haplotype = int(hp_tag)
        if not is_unmapped:
            legacy , read_info = check_read_mapping_confidence(read_id, flags, ref_id, ref_start, strand, cigar,  mapq, num_mismatches, sa_tag,MIN_ALIGNED_LENGTH, MIN_ALIGNED_RATE,
                                                 max_read_error, MAX_SEGMENTS)
            for inf in read_info:
                if inf.mapq < min_mapq:
                    for i in range(inf.ref_start//2000,inf.ref_end//2000):
                        lowmapq.append((inf.ref_id,i))                   
                else:
                    all_reads.append((inf.read_id, inf.ref_id, inf.ref_start,inf.read_end,inf.mismatch_rate))
            if legacy and inf.mapq > min_mapq:
                new_segment, ins_list = get_segment(read_id, ref_id, ref_start, strand, cigar, haplotype, mapq, genome_id, sv_size)
                for new_seg in new_segment:
                    if new_seg.read_end - new_seg.read_start >= MIN_SEGMENT_LENGTH:
                        alignments.append(new_seg)
                for ins in ins_list:
                    inslist.append(ins)
    return lowmapq,all_reads, alignments , inslist 



def get_all_reads_parallel(bam_file, num_threads, aln_dump_stream, ref_lengths, max_read_error,
                           min_mapq, genome_id,sv_size,write_reads,all_reads_bisect,lowmapq_reg,all_reads_stat):
    CHUNK_SIZE = 10000000
    print('bam_nre')
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
            
    tasks = [(bam_file, region, max_read_error, min_mapq, genome_id,sv_size) for region in fetch_list]

    parsing_results = None
    thread_pool = Pool(num_threads)
    parsing_results = thread_pool.starmap(get_split_reads, tasks)
    all_reads=[]

    segments_by_read = defaultdict(list)
    inslist = []
    lowmapq_reads = defaultdict(int)
    for lowmapq,read_stat,alignments,ins_list in parsing_results:
        for read1 in lowmapq:
            lowmapq_reads[read1]+=1
        for aln in alignments:
            segments_by_read[aln.read_id].append(aln)
        for ins in ins_list:
            inslist.append(ins)
        for allr in read_stat:
            all_reads_stat.append(allr)
            
    lowmapq_reads = dict(sorted(lowmapq_reads.items()))

    for key, value in lowmapq_reads.items():
        if value > min_lowmapq_reads:
            if not lowmapq_reg[key[0]]:
                lowmapq_reg[key[0]].append([key[1]*10000])
                lowmapq_reg[key[0]].append([(key[1]+1)*10000])
            elif (key[1]*10000)==lowmapq_reg[key[0]][1][-1]:
                lowmapq_reg[key[0]][1][-1]=(key[1]+1)*10000
            else:
                lowmapq_reg[key[0]][0].append(key[1]*10000)
                lowmapq_reg[key[0]][1].append((key[1]+1)*10000)
    for read in segments_by_read:
        segments = segments_by_read[read]
        segments.sort(key=lambda s: s.read_start)
        dedup_segments = []
        for seg in segments:
            if not dedup_segments or dedup_segments[-1].read_start != seg.read_start:
                dedup_segments.append(seg)
                if not all_reads_bisect[seg.ref_id]:
                    all_reads_bisect[seg.ref_id].append([seg.ref_start])
                    all_reads_bisect[seg.ref_id].append([seg.ref_end])
                    all_reads_bisect[seg.ref_id].append([seg.read_id])
                    all_reads_bisect[seg.ref_id].append([(seg.haplotype,seg.genome_id)])
                else:    
                    all_reads_bisect[seg.ref_id][0].append(seg.ref_start)
                    all_reads_bisect[seg.ref_id][1].append(seg.ref_end)
                    all_reads_bisect[seg.ref_id][2].append(seg.read_id)
                    all_reads_bisect[seg.ref_id][3].append((seg.haplotype,seg.genome_id))
        if len(dedup_segments)>1:
            all_reads.append(dedup_segments)
        if write_reads:
            aln_dump_stream.write(str(read) + "\n")
            for seg in dedup_segments:
                aln_dump_stream.write(str(seg) + "\n")
            aln_dump_stream.write("\n")
    
    return lowmapq_reg,all_reads_stat,all_reads,all_reads_bisect , inslist






    

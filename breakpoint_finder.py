#!/usr/bin/env python3

"""
This script takes (phased) bam file as input, and outputs coordinates
of tentative 2-breaks, in particular inversion coordinates
"""

import sys
import re
import shutil
import numpy as np
import math
from copy import copy
from collections import namedtuple, defaultdict, Counter
import pysam
from multiprocessing import Pool
import random
import os
import subprocess
import bisect
import logging


from resolve_vntr import resolve_vntr
from bam_processing import filter_all_reads, get_allsegments


logger = logging.getLogger()


MAX_LOWMAPQ_READS = 10
MIN_SEGMENT_LENGTH = 100
MIN_SEGMENT_OVERLAP = 100
MAX_SEGMENT_OVERLAP = 500
MAX_CONNECTION= 1000
MAX_UNALIGNED_LEN = 500
COV_WINDOW = 500

class ReadConnection(object):
    __slots__ = ("ref_id_1", "pos_1", "sign_1", "ref_id_2", "pos_2", "sign_2",
                 "haplotype_1", "haplotype_2", "read_id", "genome_id")
    def __init__(self, ref_id_1, pos_1, sign_1, ref_id_2, pos_2, sign_2, haplotype_1, haplotype_2, read_id, genome_id):
        self.ref_id_1 = ref_id_1
        self.ref_id_2 = ref_id_2
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.sign_1 = sign_1
        self.sign_2 = sign_2
        self.haplotype_1 = haplotype_1
        self.haplotype_2 = haplotype_2
        self.read_id = read_id
        self.genome_id = genome_id
    def signed_coord_1(self):
        return self.sign_1 * self.pos_1
    def signed_coord_2(self):
        return self.sign_2 * self.pos_2


class Breakpoint(object):
    __slots__ = ("ref_id", "position","dir_1", "spanning_reads", "connections", 
                 "read_ids", "pos2", "is_insertion", "insertion_size")
    def __init__(self, ref_id, ref_position, dir_1):
        self.ref_id = ref_id
        self.position = ref_position
        self.dir_1 = dir_1
        self.spanning_reads = defaultdict(int)
        self.connections = defaultdict(list)
        self.read_ids=[]
        self.pos2 = []
        self.is_insertion = False
        self.insertion_size = 0

    def fancy_name(self):
        if not self.is_insertion:
            return self.unique_name()
        else:
            return f"INS:{self.insertion_size}"

    def unique_name(self):
        if not self.is_insertion:
            sign = '-' if self.dir_1 == -1 else '+'
            return f"{sign}{self.ref_id}:{self.position}"
        else:
            return f"INS:{self.ref_id}:{self.position}"


class DoubleBreak(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "genome_id","haplotype_1",'haplotype_2',
                 "supp",'supp_read_ids','length','genotype','edgestyle')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, genome_id, haplotype_1, haplotype_2, 
                 supp, supp_read_ids, length, genotype, edgestyle):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
        self.direction_2 = direction_2
        self.genome_id = genome_id
        self.haplotype_1 = haplotype_1
        self.haplotype_2 = haplotype_2
        self.supp = supp
        self.supp_read_ids = supp_read_ids
        self.length = length
        self.genotype = genotype
        self.edgestyle = edgestyle
    #def directional_coord_1(self):
    #    return self.direction_1 * self.bp_1.position
    #def directional_coord_2(self):
    #    return self.direction_2 * self.bp_2.position
    def to_string(self):
        strand_1 = "+" if self.direction_1 > 0 else "-"
        strand_2 = "+" if self.direction_2 > 0 else "-"
        label_1 = "{0}{1}:{2}".format(strand_1, self.bp_1.ref_id, self.bp_1.position)
        label_2 = "{0}{1}:{2}".format(strand_2, self.bp_2.ref_id, self.bp_2.position)
        if label_2[1:] < label_1[1:]:
            label_1, label_2 = label_2, label_1
        bp_name = label_1 + "|" + label_2
        return bp_name
    

class GenomicSegment(object):
    __slots__ = "genome_id","haplotype", "ref_id", "dir1", "pos1", "dir2" , "pos2", "coverage", "length_bp"
    def __init__(self, genome_id, haplotype, ref_id, pos1, pos2, coverage, length_bp):
        self.genome_id = genome_id
        self.haplotype = haplotype
        self.ref_id = ref_id
        self.dir1 = '-'
        self.pos1 = pos1
        self.dir2 = '+'
        self.pos2 = pos2
        self.coverage = coverage
        self.length_bp = length_bp


def count_spanning_reads(read_segments, position, read_ids):
    read_segments_start = read_segments[0]
    read_segments_end = read_segments[1]
    strt = read_segments_start[1][:bisect.bisect_left(read_segments_start[0],position)]
    end = read_segments_end[1][bisect.bisect_left(read_segments_end[0],position):]
    idx = list(set.intersection(set(strt), set(end)))
    idx_filt = [i for i in idx if read_segments[4][read_segments[5].index(i)] not in read_ids]
    count_all = Counter([read_segments[3][read_segments[5].index(i)] for i in idx_filt])
    return count_all
    

def get_breakpoints(split_reads, thread_pool, ref_lengths, lowmapq_reg, args):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """
    clust_len = args.bp_cluster_size
    min_reads = args.bp_min_support
    min_ref_flank = args.min_ref_flank 
    min_mapq = args.min_mapping_quality
    MAX_SEGMENT_DIST= 500
    
    seq_breakpoints = defaultdict(list)
    def _signed_breakpoint(seg, direction):
        ref_bp, sign = None, None
        if direction == "right":
            ref_bp = seg.ref_end if seg.strand == "+" else seg.ref_start
            sign = 1 if seg.strand == "+" else -1
        elif direction == "left":
            ref_bp = seg.ref_start if seg.strand == "+" else seg.ref_end
            sign = -1 if seg.strand == "+" else 1
        return ref_bp, sign
    def _add_double(seg_1, seg_2):
        ref_bp_1, sign_1 = _signed_breakpoint(s1, "right")
        ref_bp_2, sign_2 = _signed_breakpoint(s2, "left")
        if ref_bp_1 > ref_bp_2:
            seq_breakpoints[s2.ref_id].append(ReadConnection(s2.ref_id, ref_bp_2, sign_2, s1.ref_id, ref_bp_1, sign_1,
                                                             s2.haplotype, s1.haplotype, s1.read_id, s1.genome_id))
        else:
            seq_breakpoints[s1.ref_id].append(ReadConnection(s1.ref_id, ref_bp_1, sign_1, s2.ref_id, ref_bp_2, sign_2,
                                                             s1.haplotype, s2.haplotype, s1.read_id, s1.genome_id))
            
    for read_segments in split_reads:
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            if s1.mapq >= min_mapq and s2.mapq >= min_mapq and abs(s2.read_start - s1.read_end) < MAX_SEGMENT_DIST:
                _add_double(s1, s2)
    #bp_clusters = defaultdict(list)
    bp_list = []
    for seq, bp_pos in seq_breakpoints.items():
        clusters = []
        cur_cluster = []
        bp_pos.sort(key=lambda bp: bp.pos_1)
        lowmapq = lowmapq_reg[seq]
        for rc in bp_pos:
            if cur_cluster and rc.pos_1 - cur_cluster[-1].pos_1 > clust_len: ### Add direction here
                clusters.append(cur_cluster)
                cur_cluster = [rc]
            else:
                cur_cluster.append(rc)
        if cur_cluster:
            clusters.append(cur_cluster)
        for cl in clusters:
            unique_reads = set()
            read_ids = []
            connections = defaultdict(list)
            for x in cl:
                unique_reads.add((x.read_id, (x.genome_id,x.haplotype_1)))
                read_ids.append(x.read_id)
                connections[x.ref_id_2].append(x)
            by_genome_id = defaultdict(int)
            for read in unique_reads:
                by_genome_id[read[1]] += 1
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.pos_1 for x in cl]))
                low_mapq_pass = low_mapq_check(lowmapq,position)
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and low_mapq_pass:
                    bp = Breakpoint(seq, position, x.sign_1)
                    bp.connections = connections
                    bp.read_ids = read_ids
                    bp_list.append(bp)
                    #bp_clusters[seq].append(bp)
    #bp_list = []           
    #for key, bp_cluster in bp_clusters.items():
    #    read_segments = allreads_pos[key]
    #    for bp in bp_cluster:
    #        count_all =count_spanning_reads(read_segments, bp.position, bp.read_ids) 
    #        for gen_id,counts in count_all.items():
    #            bp.spanning_reads[gen_id]=counts
    #            bp_list.append(bp)
    tasks = [(bp_1, clust_len, min_reads, lowmapq_reg,ref_lengths,min_ref_flank) for bp_1 in bp_list]
    double_break = thread_pool.starmap(get_left_break, tasks)
    double_breaks=[]
    for db in double_break:
        if db:
            double_breaks += db
    return double_breaks
    

def get_left_break(bp_1, clust_len, min_reads, lowmapq_reg,ref_lengths,min_ref_flank):
    MAX_BP_CONNECTIONS = 33 ## if segment is connected to more than two bp
    left_break = bp_1.connections
    n_conn  = 0
    clusters = []
    double_break=[]
    for seq, conn in left_break.items():
        lowmapq = lowmapq_reg[seq]
        cur_cluster = []
        conn.sort(key=lambda cn:cn.pos_2)
        for rc in conn:
            if cur_cluster and rc.pos_2 - cur_cluster[-1].pos_2 > clust_len: 
                if len(cur_cluster) > min_reads:
                    clusters.append(cur_cluster)
                n_conn +=1
                cur_cluster = [rc]
            else:
                cur_cluster.append(rc)
        if len(cur_cluster) >= min_reads:
            clusters.append(cur_cluster)
            n_conn +=1
        if n_conn > MAX_BP_CONNECTIONS or not clusters:
            continue
        for cl in clusters:
            unique_reads = defaultdict(set)
            for x in cl:
                unique_reads[(x.genome_id,x.haplotype_1,x.haplotype_2)].add(x.read_id)
            by_genome_id = defaultdict(int)
            happ_support_1 = defaultdict(list)
            happ_support_2 = defaultdict(list)
            read_ids=[]
            for key, values in unique_reads.items():
                by_genome_id[key] = len(values)
                happ_support_1[key[0]].append(key[1])
                happ_support_2[key[0]].append(key[2])
                read_ids.append(values)
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.pos_2 for x in cl]))
                low_mapq_pass = low_mapq_check(lowmapq,position)
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and low_mapq_pass:
                    bp_2 = Breakpoint(seq, position, x.sign_2)
                    #count_all =count_spanning_reads(read_segments, bp_2.position, read_ids)
                    #for gen_id,counts in count_all.items():
                    #    bp_2.spanning_reads[gen_id]=counts
                    for keys , support_reads in unique_reads.items():
                        genome_id = keys[0]
                        haplotype_1 = keys[1]
                        haplotype_2 = keys[2]
                        if sum(happ_support_1[genome_id]) == 3 or sum(happ_support_2[genome_id]) == 3:
                            genotype = 'hom'
                        else:
                            genotype = 'het'
                        supp= len(support_reads)
                        if bp_1.ref_id == bp_2.ref_id:
                            length_bp = abs(bp_1.position - bp_2.position)
                        else:
                            length_bp = 0
                        double_break.append(DoubleBreak(bp_1, bp_1.dir_1, bp_2, bp_2.dir_1,genome_id,haplotype_1, haplotype_2,supp,support_reads, length_bp, genotype , 'dashed'))
        return double_break

### Make it optional
def support_read_filter(bp_cluster, read_segments, by_genome_id, true_bp_threshold):
    LOW_SUPP_THR = 3
    HIGH_SUPP_THR= 0.75
    if bp_cluster.dir_1 == 1:
        strt=bisect.bisect_left(read_segments[1],bp_cluster.position)
        end=bisect.bisect_left(read_segments[1],bp_cluster.position+1)
    else:
        strt=bisect.bisect_left(read_segments[0],bp_cluster.position)
        end=bisect.bisect_left(read_segments[0],bp_cluster.position+1)
    count_all = Counter([read_segments[3][i] for i in range(strt,end)])
    for key,value in count_all.items():
        supp = by_genome_id[key]
        if supp<11 and value - supp >LOW_SUPP_THR:
            return False
        elif value/supp > HIGH_SUPP_THR:
            return False
    return True


def low_mapq_check(lowmapq, position):
    if lowmapq and (bisect.bisect_left(lowmapq[0], position) == 0 or lowmapq[1][bisect.bisect_left(lowmapq[0], position) - 1] < position):
        return True
        
      
def extract_insertions(ins_list_all, lowmapq_reg, clust_len, min_ref_flank, ref_lengths, min_reads):
    NUM_HAPLOTYPES = 3
    ins_list = defaultdict(list)
    for ins in ins_list_all:
        ins_list[ins.ref_id].append(ins)
   
    #t = 0
    ins_clusters = []
    for seq, ins_pos in ins_list.items():
        clusters = []
        cur_cluster = []
        ins_pos.sort(key=lambda x:x.ref_end)
        lowmapq = lowmapq_reg[seq]
        #read_segments = allreads_pos[seq]
        
        for rc in ins_pos:
            if cur_cluster and rc.ref_end - cur_cluster[-1].ref_end > clust_len:
                cur_cluster.sort(key=lambda x:x.segment_length)
                cl_ins = []
                for cl1 in cur_cluster:
                    if cl_ins and abs(cl1.segment_length - cl_ins[-1].segment_length) > clust_len:
                        clusters.append(cl_ins)
                        cl_ins = [cl1]
                    else:
                        cl_ins.append(cl1)
                if cl_ins:
                    clusters.append(cl_ins)
                cur_cluster = [rc]
            else:
                cur_cluster.append(rc)
        if cur_cluster:
            clusters.append(cur_cluster)
            
        for cl in clusters:
            unique_reads = defaultdict(set)
            for x in cl:
                unique_reads[(x.genome_id,x.haplotype)].add(x.read_id)
            by_genome_id = defaultdict(int)
            happ_support_1 = defaultdict(list)
            read_ids=[]
            for key, values in unique_reads.items():
                by_genome_id[key] = len(values)
                happ_support_1[key[0]].append(key[1])
                read_ids.append(values)
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.ref_end for x in cl]))
                ins_length = int(np.median([x.segment_length for x in cl]))
                low_mapq_pass = low_mapq_check(lowmapq,position)
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and low_mapq_pass:
                    for key, value in unique_reads.items():
                        bp_1 = Breakpoint(seq, position,-1)
                        bp_2 = Breakpoint(seq, position,1)
                        bp_1.read_ids = value
                        bp_2.read_ids = value

                        bp_3 = Breakpoint(seq, position, 1)
                        bp_3.is_insertion = True
                        bp_3.insertion_size = ins_length

                        genome_id = key[0]
                        if sum(happ_support_1[genome_id]) == NUM_HAPLOTYPES:
                            genotype = 'hom'
                        else:
                            genotype = 'het'
                        #count_all =count_spanning_reads(read_segments, position, read_ids)
                        #for gen_id,counts in count_all.items():
                        #    bp_1.spanning_reads[gen_id]=counts
                        #    bp_2.spanning_reads[gen_id]=counts
                        #    bp_3.spanning_reads[gen_id]=counts
                        db_1 = DoubleBreak(bp_1, -1, bp_3, 1, genome_id, key[1], key[1], len(value), value, ins_length, genotype, 'dashed')
                        db_2 = DoubleBreak(bp_3, 1, bp_2, 1, genome_id, key[1], key[1], len(value), value, ins_length, genotype, 'dashed')
                        ins_clusters.append(db_1)
                        ins_clusters.append(db_2)
    return(ins_clusters)


def compute_bp_coverage(double_breaks,coverage_histograms):
    for db in double_breaks:
        if not db.bp_1.is_insertion:
            db.bp_1.spanning_reads[(db.genome_id, db.haplotype_1)] = bp_coverage(db.bp_1, db.genome_id, db.haplotype_1, coverage_histograms)
        if not db.bp_2.is_insertion:
            db.bp_2.spanning_reads[(db.genome_id, db.haplotype_2)] = bp_coverage(db.bp_2, db.genome_id, db.haplotype_2, coverage_histograms)
        

def bp_coverage(bp, genome_id, haplotype, coverage_histograms):
    hist_start = bp.position // COV_WINDOW
    cov_bp = coverage_histograms[(genome_id, haplotype, bp.ref_id)][hist_start]
    if not cov_bp:
        return 0
    return cov_bp
    

def get_phasingblocks(hb_vcf):
    MIN_BLOCK_LEN = 10000
    MIN_SNP = 10

    vcf = pysam.VariantFile(hb_vcf)
    haplotype_blocks = defaultdict(list)
    endpoint_list = defaultdict(list)
    switch_points = defaultdict(list)

    for var in vcf:
        if 'PS' in var.samples.items()[0][1].items()[-1]:
            haplotype_blocks[(var.chrom, var.samples.items()[0][1]['PS'])].append(var.pos)

    for (chr_id, block_name), coords in haplotype_blocks.items():
        if max(coords) - min(coords) > MIN_BLOCK_LEN and len(coords) >= MIN_SNP:
            endpoint_list[chr_id].append(min(coords))
            endpoint_list[chr_id].append(max(coords))
            #print(len(coords), min(coords), max(coords))

    for chr_id, values in endpoint_list.items():
        values.sort()
        switch_points[chr_id] = [(a + b) // 2 for a, b in zip(values[:-1], values[1:])]

    return switch_points


def segment_coverage(histograms, genome_id, ref_id, ref_start, ref_end, haplotype):
    hist_start = ref_start // COV_WINDOW
    hist_end = ref_end // COV_WINDOW
    cov_list = histograms[(genome_id, haplotype, ref_id)][hist_start : hist_end + 1]

    if not cov_list:
        return 0
    return int(np.median(cov_list))


def get_segments_coverage(segments, coverage_histograms):
    genomic_segments = []
    for (genome_id, seg_ref, seg_start, seg_end, seg_hp) in segments:
        coverage = segment_coverage(coverage_histograms, genome_id, seg_ref, seg_start, seg_end, seg_hp)
        genomic_segments.append(GenomicSegment(genome_id, seg_hp, seg_ref, seg_start, seg_end,
                                               coverage, seg_end - seg_start))

    return genomic_segments


def get_genomic_segments(double_breaks, coverage_histograms, thread_pool, hb_vcf):
    switch_points = defaultdict(list)
    if hb_vcf:
        switch_points = get_phasingblocks(hb_vcf)

    single_bp = defaultdict(list)
    genomic_segments=[]
    segments = []
    for double_bp in double_breaks:
        if not double_bp.bp_1.is_insertion:
            single_bp[(double_bp.genome_id, double_bp.haplotype_1, double_bp.bp_1.ref_id)].append(double_bp.bp_1.position)
        if not double_bp.bp_2.is_insertion:
            single_bp[(double_bp.genome_id, double_bp.haplotype_2, double_bp.bp_2.ref_id)].append(double_bp.bp_2.position)          

    for (genome_name, haplotype_name, ref_name), s_bp in single_bp.items():
        s_bp1 = s_bp + switch_points[ref_name]
        s_bp1 = list(set(s_bp1))
        s_bp1.sort()
        for seg_start, seg_end in zip(s_bp1[:-1], s_bp1[1:]):
            segments.append((genome_name, ref_name, seg_start, seg_end, haplotype_name))

    genomic_segments = get_segments_coverage(segments, coverage_histograms)
    return genomic_segments, switch_points


def extract_lowmapq_regions(segments_by_read , min_mapq):
    MAX_LOWMAPQ_READS = 10
    MIN_LOWMAPQ_RATIO = 0.5
    WIN_SIZE = 2000
        
    lowmapq_reads = defaultdict(int)
    highmapq_reads = defaultdict(int)
    for segments in segments_by_read.values():
        segments.sort(key=lambda s: s.read_start)
        for seg in segments:
            if seg.mapq < min_mapq and not seg.is_insertion:
                for i in range(seg.ref_start//WIN_SIZE,seg.ref_end//WIN_SIZE): ###Insertion
                    lowmapq_reads[(seg.ref_id,i)]+=1
            elif not seg.is_insertion:
                for i in range(seg.ref_start//WIN_SIZE,seg.ref_end//WIN_SIZE): ###Insertion
                    highmapq_reads[(seg.ref_id,i)]+=1
                    
    lowmapq_reads = dict(sorted(lowmapq_reads.items()))
    lowmapq_reg = defaultdict(list)
    for key, value in lowmapq_reads.items():
        if value > highmapq_reads[key]*MIN_LOWMAPQ_RATIO:
            if not lowmapq_reg[key[0]]:
                lowmapq_reg[key[0]].append([key[1]*WIN_SIZE])
                lowmapq_reg[key[0]].append([(key[1]+1)*WIN_SIZE])
            elif (key[1]*WIN_SIZE)==lowmapq_reg[key[0]][1][-1]:
                lowmapq_reg[key[0]][1][-1]=(key[1]+1)*WIN_SIZE
            else:
                lowmapq_reg[key[0]][0].append(key[1]*WIN_SIZE)
                lowmapq_reg[key[0]][1].append((key[1]+1)*WIN_SIZE)
    return lowmapq_reg



def get_insertionreads(segments_by_read_filtered):
    ins_list_all = []
    for read in segments_by_read_filtered:
        for seg in read:
            if seg.is_insertion:
                ins_list_all.append(seg)
    return ins_list_all
    

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


def write_alignments(allsegments, outpath):
    aln_dump_stream = open(outpath, "w")
    for read in allsegments:
        aln_dump_stream.write(str(read) + "\n")
    aln_dump_stream.write("\n")
        
def output_breaks(double_breaks, genome_tags, phasing, out_stream):
    loc = defaultdict(int)
    t = 0
    header = '#BP_pos1:BP_pos2,'
    def_array = []
    if phasing:
        hp_list = [0,1,2]
        for tag in genome_tags:
            for k in hp_list:
                def_array.append((0,0,0))
                loc[(tag,k)] = t
                header += '_'.join([tag, str(k),'support,'])
                header += '_'.join([tag, str(k),'spanning_1,'])
                header += '_'.join([tag, str(k),'spanning_2,'])
                t += 1
    else:
        for tag in genome_tags:
            def_array.append((0,0,0))
            loc[(tag,0)] = t
            header += '_'.join([tag, 'support,'])
            header += '_'.join([tag, 'spanning_1,'])
            header += '_'.join([tag, 'spanning_2,'])
            t += 1
    summary_csv = defaultdict(list)
    for br in double_breaks:
        if not br.bp_1.is_insertion:
            if not summary_csv[br.to_string()]:
                summary_csv[br.to_string()] = def_array[:]
                idd=(br.genome_id, br.haplotype_1)
                summary_csv[br.to_string()][loc[idd]] = (br.supp, br.bp_1.spanning_reads[idd], br.bp_2.spanning_reads[idd])            
            else:
                idd=(br.genome_id, br.haplotype_1)
                summary_csv[br.to_string()][loc[idd]] = (br.supp, br.bp_1.spanning_reads[idd], br.bp_2.spanning_reads[idd])
    out_stream.write(header + "\n")
    for key,values in summary_csv.items():
        bp_array=[]
        for k in values:
            bp_array.append(str(k[0]))
            bp_array.append(str(k[1]))
            bp_array.append(str(k[2]))
        bp_to_write = ','.join([key, ','.join(bp_array)])
        out_stream.write(bp_to_write)
        out_stream.write("\n")

        
def call_breakpoints(segments_by_read, thread_pool, ref_lengths, coverage_histograms, args):
    logger.info('Detecting low mapping quality regions')
    lowmapq_reg = extract_lowmapq_regions(segments_by_read, args.min_mapping_quality)
    logger.info('Filtering reads')
    if args.vntr_file:
        segments_by_read = resolve_vntr(segments_by_read, args.vntr_file, args.sv_size)
    segments_by_read_filtered = filter_all_reads(segments_by_read,args.min_mapping_quality, args.max_read_error)
    
    allsegments = get_allsegments(segments_by_read_filtered)
    if args.write_alignments:
        outpath_alignments = os.path.join(args.out_dir, "read_alignments")
        write_alignments(allsegments, outpath_alignments)
    #allreads_pos = all_reads_position(allsegments)
    
    split_reads = get_splitreads(segments_by_read_filtered)
    logger.info('Resolving overlaps')
    split_reads = resolve_overlaps(split_reads,  args.sv_size)
    
    ins_list_all = get_insertionreads(segments_by_read_filtered)
   
    logger.info('Starting breakpoint detection')
    double_breaks = get_breakpoints(split_reads,thread_pool, ref_lengths, lowmapq_reg, args)
    double_breaks.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position, b.direction_1))
    
    logger.info('Clustering unmapped insertions')
    ins_clusters = extract_insertions(ins_list_all, lowmapq_reg,  args.bp_cluster_size, args.min_ref_flank, ref_lengths,args.bp_min_support)
    double_breaks +=  ins_clusters 
            
    compute_bp_coverage(double_breaks,coverage_histograms)
    return double_breaks 


    
    

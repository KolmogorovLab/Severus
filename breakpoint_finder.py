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


MAX_LOWMAPQ_READS = 10
MIN_SEGMENT_LENGTH = 100
MIN_SEGMENT_OVERLAP = 100
MAX_SEGMENT_OVERLAP = 500
MAX_CONNECTION= 1000

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
    __slots__ = "ref_id", "position","dir_1", "spanning_reads", "connections" , "read_ids", "pos2"
    def __init__(self, ref_id, ref_position,dir_1):
        self.ref_id = ref_id
        self.position = ref_position
        self.dir_1 = dir_1
        self.spanning_reads = defaultdict(int)
        self.connections =defaultdict(list)
        self.read_ids=[]
        self.pos2 = []


class DoubleBreak(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "genome_id","haplotype1",'haplotype2',
                 "supp",'supp_read_ids','length','genotype','edgestyle')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, genome_id, haplotype1, haplotype2, 
                 supp, supp_read_ids, length, genotype, edgestyle):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
        self.direction_2 = direction_2
        self.genome_id = genome_id
        self.haplotype1 = haplotype1
        self.haplotype2 = haplotype2
        self.supp = supp
        self.supp_read_ids = supp_read_ids
        self.length = length
        self.genotype = genotype
        self.edgestyle = edgestyle
    def directional_coord_1(self):
        return self.direction_1 * self.bp_1.position
    def directional_coord_2(self):
        return self.direction_2 * self.bp_2.position
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

def read_vntr_file(vntr_file):
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

def resolve_vntrs(seq_breakpoints,vntr_list):
    for seq, bp_pos in seq_breakpoints.items():
        tr_reg = vntr_list[seq]
        for rc in bp_pos:
            if rc.ref_id_1 == rc.ref_id_2:
                strt = bisect.bisect_left(tr_reg[0],rc.pos_1)
                end = bisect.bisect_left(tr_reg[1],rc.pos_2)
                if strt - end == 1:
                    length = abs(rc.pos_1 - rc.pos_2)
                    if rc.pos_1 < rc.pos_2:
                        rc.pos_1 = tr_reg[0][end]
                        rc.pos_2 = rc.pos_1 + length
                    else:
                        rc.pos_2 = tr_reg[0][end]
                        rc.pos_1 = rc.pos_2 + length
    return seq_breakpoints

def get_breakpoints(allreads_pos, split_reads,vntr_list, thread_pool,clust_len, max_unaligned_len,
                    min_reads, min_ref_flank, ref_lengths, min_mapq,single_bp,lowmapq_reg,min_sv_size):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """
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
        seq_breakpoints[s1.ref_id].append(ReadConnection(s1.ref_id, ref_bp_1, sign_1, s2.ref_id, ref_bp_2, sign_2,
                                                         s1.haplotype, s2.haplotype, s1.read_id, s1.genome_id))
    def _add_single(seg, direction):
        ref_bp, sign = _signed_breakpoint(seg, direction)
        seq_breakpoints[seg.ref_id].append(ReadConnection(seg.ref_id, ref_bp, sign, None, None, None,
                                                          seg.haplotype, None, seg.read_id, seg.genome_id))
    
    for read_segments in split_reads:
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            if s1.mapq >= min_mapq and s2.mapq >= min_mapq:
                _add_double(s1, s2)
            elif single_bp:
                if s1.mapq >= min_mapq:
                    _add_single(s1, "right")
                if s2.mapq >= min_mapq:
                    _add_single(s2, "left")
    if vntr_list:
        seq_breakpoints = resolve_vntrs(seq_breakpoints,vntr_list)
                    
    bp_clusters = defaultdict(list)
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
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and lowmapq and lowmapq[1][bisect.bisect_left(lowmapq[0],position)-1]<position:
                    bp_cluster = Breakpoint(seq, position, x.sign_1)
                    bp_cluster.connections = connections
                    bp_cluster.read_ids = read_ids
                    bp_clusters[seq].append(bp_cluster)
    bp_list = []           
    for key, bp_cluster in bp_clusters.items():
        read_segments = allreads_pos[key]
        read_segments[0], read_segments[1], read_segments[2], read_segments[3] = zip(*sorted(zip(read_segments[0],read_segments[1],read_segments[2],read_segments[3])))
        for bp in bp_cluster:
            strt=bisect.bisect_left(read_segments[1],bp.position)
            end = bisect.bisect_left(read_segments[0],bp.position)
            count_all = Counter([read_segments[3][i] for i in range(strt,end) if read_segments[1][i] > bp.position]) 
            for gen_id,counts in count_all.items():
                bp.spanning_reads[gen_id]=counts
                bp_list.append(bp)
    
    tasks = [(bp_1, clust_len, min_reads, lowmapq_reg,ref_lengths,min_ref_flank,allreads_pos) for bp_1 in bp_list]
    double_break = thread_pool.starmap(get_left_break, tasks)
         
    double_breaks=[]
    for db in double_break:
        if db:
            double_breaks += db
    return double_breaks
    

def get_left_break(bp_1, clust_len, min_reads, lowmapq_reg,ref_lengths,min_ref_flank,allreads_pos):
    MAX_BP_CONNECTIONS = 3 ## if segment is connected to more than two bp
    left_break = bp_1.connections
    n_conn  = 0
    clusters = []
    double_break=[]
    for seq, conn in left_break.items():
        lowmapq = lowmapq_reg[seq]
        read_segments = allreads_pos[seq]
        read_segments[0], read_segments[1], read_segments[2], read_segments[3] = zip(*sorted(zip(read_segments[0],read_segments[1],read_segments[2],read_segments[3])))
        
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
        if len(cur_cluster) > min_reads:
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
            for key, values in unique_reads.items():
                by_genome_id[key] = len(values)
                happ_support_1[key[0]].append(key[1])
                happ_support_2[key[0]].append(key[2])
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.pos_2 for x in cl]))
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and lowmapq and lowmapq[1][bisect.bisect_left(lowmapq[0],position)-1]<position:
                    bp_2 = Breakpoint(seq, position, x.sign_2)
                    strt=bisect.bisect_left(read_segments[1],bp_2.position)
                    end = bisect.bisect_left(read_segments[0],bp_2.position)
                    count_all = Counter([read_segments[3][i] for i in range(strt,end) if read_segments[1][i] > bp_2.position]) 
                    for gen_id,counts in count_all.items():
                        bp_2.spanning_reads[gen_id]=counts
                    for keys , support_reads in unique_reads.items():
                        genome_id = keys[0]
                        haplotype1 = keys[1]
                        haplotype2 = keys[2]
                        if sum(happ_support_1[genome_id]) == 3 or sum(happ_support_2[genome_id]) == 3:
                            genotype = 'hom'
                        else:
                            genotype = 'het'
                        supp= len(support_reads)
                        if bp_1.ref_id == bp_2.ref_id:
                            length_bp = abs(bp_1.position - bp_2.position)
                        else:
                            length_bp = 0
                        double_break.append(DoubleBreak(bp_1, bp_1.dir_1, bp_2, bp_2.dir_1,genome_id,haplotype1, haplotype2,supp,support_reads, length_bp, genotype , 'dashed'))
        return double_break


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


def resolve_vntr_ins(ins_list,vntr_list):
    for seq, reads in ins_list.items():
        tr_reg = vntr_list[seq]
        for rc in reads:
            strt = bisect.bisect_right(tr_reg[0],rc.ref_end)
            end = bisect.bisect_left(tr_reg[1],rc.ref_end)
            if strt - end == 1:
                rc.ref_end = tr_reg[0][end]
    return ins_list

      
def extract_insertions(ins_list_all, lowmapq_reg, vntr_list, clust_len, min_ref_flank, ref_lengths, min_reads, allreads_pos):
    ins_list = defaultdict(list)
    for ins in ins_list_all:
        ins_list[ins.ref_id].append(ins)
    
    if vntr_list:
        ins_list = resolve_vntr_ins(ins_list, vntr_list)
   
    #t = 0
    ins_clusters = []
    for seq, ins_pos in ins_list.items():
        clusters = []
        cur_cluster = []
        ins_pos.sort(key=lambda x:x.ref_end)
        lowmapq = lowmapq_reg[seq]
        read_segments = allreads_pos[seq]
        read_segments[0], read_segments[1], read_segments[2], read_segments[3] = zip(*sorted(zip(read_segments[0],read_segments[1],read_segments[2],read_segments[3])))
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
                    clusters.append( cl_ins)
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
            for key, values in unique_reads.items():
                by_genome_id[key] = len(values)
                happ_support_1[key[0]].append(key[1])
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.ref_end for x in cl]))
                ins_length = int(np.median([x.segment_length for x in cl]))
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and lowmapq and lowmapq[1][bisect.bisect_left(lowmapq[0],position)-1]<position:
                    #t+=1
                    for key, value in unique_reads.items():
                        bp_1 = Breakpoint(seq, position,-1)
                        bp_2 = Breakpoint(seq, position,1)
                        bp_1.read_ids = value
                        bp_2.read_ids = value
                        bp_3 = Breakpoint('INS', ins_length,1)
                        genome_id = key[0]
                        if sum(happ_support_1[genome_id]) == 3:
                            genotype = 'hom'
                        else:
                            genotype = 'het'
                        count_all = Counter(read_segments[3][bisect.bisect_left(read_segments[1],position):bisect.bisect_left(read_segments[0],position)])
                        for gen_id,counts in count_all.items():
                            bp_1.spanning_reads[gen_id]=counts-len(value)
                            bp_2.spanning_reads[gen_id]=counts-len(value)
                            bp_3.spanning_reads[gen_id]=counts-len(value)
                        db = DoubleBreak(bp_1, -1, bp_3, 1, genome_id, key[1], key[1], len(value), value, ins_length, genotype, 'dashed')
                        db = DoubleBreak(bp_3, 1, bp_2, 1, genome_id, key[1], key[1], len(value), value, ins_length, genotype, 'dashed')
                        ins_clusters.append(db)
    return(ins_clusters)


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
    COV_WINDOW = 500
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
        if not double_bp.bp_1.ref_id == 'INS' and not double_bp.bp_2.ref_id == 'INS':
            single_bp[(double_bp.genome_id, double_bp.haplotype1, double_bp.bp_1.ref_id)].append(double_bp.bp_1.position)
            single_bp[(double_bp.genome_id, double_bp.haplotype2, double_bp.bp_2.ref_id)].append(double_bp.bp_2.position)

    for (genome_name, haplotype_name, ref_name), s_bp in single_bp.items():
        s_bp1 = s_bp + switch_points[ref_name]
        s_bp1 = list(set(s_bp1))
        s_bp1.sort()
        for seg_start, seg_end in zip(s_bp1[:-1], s_bp1[1:]):
            segments.append((genome_name, ref_name, seg_start, seg_end, haplotype_name))

    genomic_segments = get_segments_coverage(segments, coverage_histograms)
    return genomic_segments, switch_points
    
def output_breaks(double_breaks, genome_tags,phasing, out_stream):
    loc = defaultdict(int)
    t= 0
    header = '#BP_pos1:BP_pos2,'
    def_array = []
    hp_list = [0,1,2]
    if phasing:
        for tag in genome_tags:
            for k in hp_list:
                def_array.append((0,0,0))
                loc[(k,tag)] = t
                header += '_'.join([tag, str(k),'support,'])
                header += '_'.join([tag, str(k),'spanning_1,'])
                header += '_'.join([tag, str(k),'spanning_2,'])
                t+=1
    else:
        for tag in genome_tags:
            def_array.append((0,0,0))
            loc[(0,tag)] = t
            header += '_'.join([tag, 'support,'])
            header += '_'.join([tag, 'spanning_1,'])
            header += '_'.join([tag, 'spanning_2,'])
            t+=1
            
    summary_csv = defaultdict(list)
    for br in double_breaks:
        if not summary_csv[br.to_string()]:
            summary_csv[br.to_string()] = def_array[:]
            idd=(br.haplotype1,br.genome_id)
            summary_csv[br.to_string()][loc[idd]] = (br.supp , br.bp_1.spanning_reads[idd],br.bp_2.spanning_reads[idd])            
        else:
            idd=(br.haplotype1,br.genome_id)
            summary_csv[br.to_string()][loc[idd]] = (br.supp , br.bp_1.spanning_reads[idd],br.bp_2.spanning_reads[idd])
            
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
        
        
def call_breakpoints(split_reads,allreads_pos,ins_list_all,lowmapq_reg,vntr_file,thread_pool,clust_len,max_unaligned_len,min_reads, min_ref_flank, ref_lengths, min_mapq,single_bp,min_sv_size,max_read_error):
    vntr_list = []
    if vntr_file:
        vntr_list=read_vntr_file(vntr_file)
    print('Starting breakpoint detection')
    double_breaks = get_breakpoints(allreads_pos, split_reads,vntr_list, thread_pool,clust_len, max_unaligned_len,min_reads, min_ref_flank, ref_lengths, min_mapq,single_bp,lowmapq_reg,min_sv_size)
    print('Clustering unmapped insertions')
    ins_clusters = extract_insertions(ins_list_all, lowmapq_reg,vntr_list,clust_len,min_ref_flank,ref_lengths,min_reads,allreads_pos)
    double_breaks +=  ins_clusters             
    double_breaks.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position, b.direction_1))
    return double_breaks 


    
    

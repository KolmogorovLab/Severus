#!/usr/bin/env python3

"""
This script takes (phased) bam file as input, and outputs coordinates
of tentative 2-breaks, in particular inversion coordinates
"""


import numpy as np
from collections import defaultdict, Counter
import pysam
import os
import bisect
import logging
import networkx as nx
import copy

from severus.bam_processing import _calc_nx, extract_clipped_end, get_coverage_parallel

logger = logging.getLogger()


MAX_LOWMAPQ_READS = 10
MIN_SEGMENT_LENGTH = 100
MIN_SEGMENT_OVERLAP = 100
MAX_SEGMENT_OVERLAP = 500
MAX_CONNECTION= 1000
MAX_UNALIGNED_LEN = 500
COV_WINDOW = 500

class ReadConnection(object):
    __slots__ = ("ref_id_1", "pos_1", "pos_1_ori", "sign_1", "ref_id_2", "pos_2", "pos_2_ori", "sign_2", "haplotype_1", "haplotype_2", "read_id", 
                 "genome_id", 'bp_list', 'is_pass1', 'is_pass2', 'mapq_1', 'mapq_2', 'is_dup', 'has_ins', 'ins_pos', 'seg_len')
    def __init__(self, ref_id_1, pos_1, pos_1_ori , sign_1, ref_id_2, pos_2, pos_2_ori, sign_2, haplotype_1, 
                 haplotype_2, read_id, genome_id, is_pass1, is_pass2, mapq_1, mapq_2, is_dup, has_ins, ins_pos, seg_len):
        self.ref_id_1 = ref_id_1
        self.ref_id_2 = ref_id_2
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.pos_1_ori = pos_1_ori
        self.pos_2_ori = pos_2_ori
        self.sign_1 = sign_1
        self.sign_2 = sign_2
        self.haplotype_1 = haplotype_1
        self.haplotype_2 = haplotype_2
        self.read_id = read_id
        self.genome_id = genome_id
        self.bp_list = []
        self.is_pass1 = is_pass1
        self.is_pass2 = is_pass2
        self.mapq_1 = mapq_1
        self.mapq_2 = mapq_2
        self.is_dup = is_dup
        self.has_ins = has_ins
        self.ins_pos = ins_pos
        self.seg_len = seg_len
    def signed_coord_1(self):
        return self.sign_1 * self.pos_1
    def signed_coord_2(self):
        return self.sign_2 * self.pos_2
    def get_pos(self, bp_dir):
        return self.pos_1 if bp_dir == 'right' else self.pos_2
    def get_pos_ori(self, bp_dir):
        return self.pos_1_ori if bp_dir == 'right' else self.pos_2_ori
    def get_qual(self, bp_dir):
        return self.mapq_1 if bp_dir == 'right' else self.mapq_2
    def is_pass(self, bp_dir):
        return self.is_pass1 if bp_dir == 'right' else self.is_pass2
    def dir_1(self, bp_dir):
        return self.sign_1 if bp_dir == 'right' else self.sign_2
    
class Breakpoint(object):
    __slots__ = ("ref_id", "position","dir_1", "spanning_reads", "connections", 'prec',
                 "read_ids", "pos2", 'id', "is_insertion", "insertion_size", "qual", 'contig_id', 'loose_end_id')
    def __init__(self, ref_id, ref_position, dir_1, qual):
        self.ref_id = ref_id
        self.position = ref_position
        self.dir_1 = dir_1
        self.spanning_reads = defaultdict(int)
        self.connections = defaultdict(list)
        self.read_ids=[]
        self.pos2 = []
        self.id = 0
        self.is_insertion = False
        self.insertion_size = None
        self.qual = qual
        self.contig_id = False
        self.loose_end_id = False
        self.prec = 1

    def fancy_name(self):
        if not self.is_insertion or self.contig_id:
            return self.unique_name()
        if self.loose_end_id:
            return ''
        else:
            return f"INS:{self.insertion_size}"

    def unique_name(self):
        if self.is_insertion:
            return f"INS:{self.ref_id}:{self.position}"
        elif self.loose_end_id:
            return f"loose:{self.ref_id}:{self.position}"
        else:
            return f"{self.ref_id}:{self.position}"

    def coord_tuple(self):
        sign = '-' if self.dir_1 == -1 else '+'
        return (self.ref_id, self.position, sign)


class DoubleBreak(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "genome_id","haplotype_1",'haplotype_2',"supp",'supp_read_ids',
                 'length','genotype','edgestyle', 'is_pass', 'ins_seq', 'mut_type', 'is_dup', 'has_ins','subgraph_id', 'sv_type',
                 'DR', 'DV', 'hvaf', 'vaf', 'prec', 'phaseset_id', 'cluster_id', 'vcf_id', 'vcf_sv_type','vaf_pass', 'vcf_qual')
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
        self.is_pass = 'PASS'
        self.ins_seq = None
        self.is_dup = None
        self.mut_type = 'germline'
        self.has_ins = 0
        self.subgraph_id = []
        self.sv_type = None
        self.DR =0
        self.DV = 0
        self.hvaf = 0.0
        self.vaf = 0.0
        self.prec = 1
        self.phaseset_id = (0,0)
        self.cluster_id = 0
        self.vcf_id = None
        self.vcf_sv_type = None
        self.vaf_pass = None
        self.vcf_qual = None
        
    def to_string(self):
        strand_1 = "+" if self.direction_1 > 0 else "-"
        strand_2 = "+" if self.direction_2 > 0 else "-"
        label_1 = "{0}{1}:{2}".format(strand_1, self.bp_1.ref_id, self.bp_1.position)
        if self.bp_2.is_insertion:
            label_1 = "{0}:{1}".format(self.bp_1.ref_id, self.bp_1.position)
            label_2 = "{0}:{1}".format('INS', self.length)
        elif self.bp_2.contig_id or self.bp_2.loose_end_id:
            label_2 = "{0}:{1}".format(self.bp_2.ref_id, self.bp_2.position)
        else:
            label_2 = "{0}{1}:{2}".format(strand_2, self.bp_2.ref_id, self.bp_2.position)
            if label_2[1:] < label_1[1:]:
                label_1, label_2 = label_2, label_1
        bp_name = label_1 + "|" + label_2
        return bp_name
    
    def to_string_csv(self):
        strand_1 = "+" if self.direction_1 > 0 else "-"
        strand_2 = "+" if self.direction_2 > 0 else "-"
        label_1 = "{0}{1}:{2}".format(strand_1, self.bp_1.ref_id, self.bp_1.position)
        if self.bp_2.is_insertion:
            label_1 = "{0}:{1}".format(self.bp_1.ref_id, self.bp_1.position)
            label_2 = "{0}:{1}".format('INS', self.length)
        elif self.bp_2.contig_id or self.bp_2.loose_end_id:
            label_2 = "{0}:{1}".format(self.bp_2.ref_id, self.bp_2.position)
        else:
            label_2 = "{0}{1}:{2}".format(strand_2, self.bp_2.ref_id, self.bp_2.position)
            if label_2[1:] < label_1[1:]:
                label_1, label_2 = label_2, label_1
        bp_name = f"{label_1}\t{label_2}"
        return bp_name
        
    

class GenomicSegment(object):
    __slots__ = "genome_id","haplotype", "ref_id", 'dir1', "pos1", 'dir2', "pos2", "coverage", "length_bp", "is_insertion"
    def __init__(self, genome_id, haplotype, ref_id, pos1, pos2, coverage, length_bp):
        self.genome_id = genome_id
        self.haplotype = haplotype
        self.ref_id = ref_id
        self.dir1 = -1
        self.dir2 = 1
        self.pos1 = pos1
        self.pos2 = pos2
        self.coverage = coverage
        self.length_bp = length_bp
        self.is_insertion = False

    def left_coord_str(self):
        return f"{self.ref_id}:{self.pos1}"

    def right_coord_str(self):
        return f"{self.ref_id}:{self.pos2}"

    def left_coord_tuple(self):
        return (self.ref_id, self.pos1, self.dir1)

    def right_coord_tuple(self):
        return (self.ref_id, self.pos2, self.dir2)
    
    def full_name(self):
        if self.is_insertion:
            return f"{self.ref_id}:{self.pos1}:INS:{self.length_bp}"
        if self.length_bp:
            return f"{self.ref_id}:{self.pos1}-{self.pos2}"
        return "swp"
    def ins_label(self):
        return f"INS:{self.length_bp}"
        
    
def get_breakpoints(split_reads, ref_lengths, args):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """

    clust_len = args.bp_cluster_size
    min_reads = args.bp_min_support
    min_ref_flank = args.min_ref_flank 
    sv_size = args.min_sv_size
    MAX_SEGMENT_DIST= 1000
    
    seq_breakpoints_l = defaultdict(list)
    seq_breakpoints_r = defaultdict(list)
    
    def _signed_breakpoint(seg, direction):
        ref_bp, sign = None, None
        if direction == "right":
            ref_bp = seg.ref_end if seg.strand == 1 else seg.ref_start
            ref_bp_ori = seg.ref_end_ori if seg.strand == 1 else seg.ref_start_ori
            sign = 1 if seg.strand == 1 else -1
        elif direction == "left":
            ref_bp = seg.ref_start if seg.strand == 1 else seg.ref_end
            ref_bp_ori = seg.ref_start_ori if seg.strand == 1 else seg.ref_end_ori
            sign = -1 if seg.strand == 1 else 1
        return ref_bp, ref_bp_ori, sign
    
    def _add_double(s1, s2):
        has_ins = 0
        dist = s2.read_start - s1.read_end
        if dist > sv_size:
            has_ins = dist
        ref_bp_1, ref_bp_1_ori, sign_1 = _signed_breakpoint(s1, "right")
        ref_bp_2, ref_bp_2_ori, sign_2 = _signed_breakpoint(s2, "left")
        is_dup = False
        if ref_bp_1 > ref_bp_2:
            if sign_1 == 1 and sign_2 == -1 and s1.ref_start <= s2.ref_start <= s1.ref_end <= s2.ref_end:
                is_dup = True
            rc = ReadConnection(s2.ref_id, ref_bp_2, ref_bp_2_ori, sign_2, s1.ref_id, ref_bp_1, ref_bp_1_ori, sign_1,
                                s2.haplotype, s1.haplotype, s1.read_id, s1.genome_id, s2.is_pass, 
                                s1.is_pass, s2.mapq, s1.mapq, is_dup, has_ins, s2.read_start,(s2.segment_length, s1.segment_length))
            seq_breakpoints_r[s2.ref_id].append(rc)
            seq_breakpoints_l[s1.ref_id].append(rc)
        else:
            if sign_1 == -1 and sign_2 == 1 and s2.ref_start <= s1.ref_start <= s2.ref_end <= s1.ref_end:
                is_dup = True
            rc = ReadConnection(s1.ref_id, ref_bp_1, ref_bp_1_ori, sign_1, s2.ref_id, ref_bp_2, ref_bp_2_ori, sign_2,
                                s1.haplotype, s2.haplotype, s1.read_id, s1.genome_id, s1.is_pass, 
                                s2.is_pass, s1.mapq, s1.mapq, is_dup, has_ins, s2.read_start, (s1.segment_length, s2.segment_length))
            seq_breakpoints_r[s1.ref_id].append(rc)
            seq_breakpoints_l[s2.ref_id].append(rc)
        
    for read_segments in split_reads:
        read_segments.sort(key=lambda x:(x.align_start, x.read_start))
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            if s2.read_start - s1.read_end < MAX_SEGMENT_DIST:
                _add_double(s1, s2)
                
    all_breaks = []
    for seq, bp_pos in seq_breakpoints_r.items():
        bps = cluster_bp(seq, bp_pos, clust_len, min_ref_flank, ref_lengths, min_reads,'right')
        if bps:
            all_breaks += bps
            
    for seq, bp_pos in seq_breakpoints_l.items():  
        bps = cluster_bp(seq, bp_pos, clust_len, min_ref_flank, ref_lengths, min_reads,'left')
        if bps:
            all_breaks += bps
            
    for bp in all_breaks:
        for conn in bp.connections:
            conn.bp_list.append(bp)
            
    matched_bp = match_breaks(seq_breakpoints_r)
    
    double_breaks=[]
    for (bp_1 , bp_2), cl in matched_bp.items():
        db = get_double_breaks(bp_1, bp_2, cl, sv_size, min_reads)
        if db:
            double_breaks += db
    double_breaks = match_del(double_breaks)
    return double_breaks


def cluster_bp(seq, bp_pos, clust_len, min_ref_flank, ref_lengths, min_reads, bp_dir):
    
    clusters = []
    cur_cluster = []
    bp_list = []
    
    bp_pos.sort(key=lambda bp: (bp.dir_1(bp_dir), bp.get_pos(bp_dir)))
    for rc in bp_pos:
        if cur_cluster and rc.get_pos(bp_dir) - cur_cluster[-1].get_pos(bp_dir) > clust_len: 
            clusters.append(cur_cluster)
            cur_cluster = [rc]
        else:
            cur_cluster.append(rc)
    if cur_cluster:
        clusters.append(cur_cluster)
        
    for cl in clusters:
        unique_reads = set()
        read_ids = []
        connections =[]
        
        for x in cl:
            unique_reads.add((x.read_id, (x.genome_id,x.haplotype_1)))
            read_ids.append(x.read_id)
            connections.append(x)
            
        by_genome_id = defaultdict(int)
        for read in unique_reads:
            by_genome_id[read[1]] += 1
            
        position_arr = [x.get_pos_ori(bp_dir) for x in cl if x.is_pass(bp_dir) == 'PASS']
        qual_arr = [x.get_qual(bp_dir) for x in cl if x.is_pass(bp_dir) == 'PASS']
        if not position_arr:
            continue
        position = int(np.median(position_arr))
        prec = 1
        if max(position_arr) - min(position_arr) > 2 * clust_len:
            prec = 0
        qual = int(np.median(qual_arr))
        sign  = x.sign_1 if bp_dir == 'right' else x.sign_2
        if position >= min_ref_flank and position <= ref_lengths[seq] - min_ref_flank:
            bp = Breakpoint(seq, position, sign, qual)
            bp.connections = connections
            bp.read_ids = read_ids
            bp.prec = prec
            bp_list.append(bp)
                
    return bp_list
            
def match_breaks(seq_breakpoints_r):
    matched_bp = defaultdict(list)
    for rc_list in seq_breakpoints_r.values():
        for rc in rc_list:
            if not len(rc.bp_list) == 2:
                continue
            rc.bp_list.sort(key=lambda bp: bp.position)
            matched_bp[(rc.bp_list[0], rc.bp_list[1])].append(rc)
    return matched_bp
        
def get_double_breaks(bp_1, bp_2, cl, sv_size, min_reads):  
    unique_reads = defaultdict(set)
    unique_reads_pass = defaultdict(set)
    db_list = []
    is_dup = []
    
    for x in cl:
        unique_reads[(x.genome_id,x.haplotype_1,x.haplotype_2)].add(x.read_id)
        is_dup.append(x.is_dup)
        if x.is_pass1 == 'PASS' and x.is_pass2 == 'PASS':
            unique_reads_pass[(x.genome_id,x.haplotype_1,x.haplotype_2)].add(x.read_id)
    
    pos1 = int(np.median([c.pos_1_ori for c in cl]))
    pos2 = int(np.median([c.pos_2_ori for c in cl]))
    if not pos1 == bp_1.position:
        bp_1 = copy.copy(bp_1)
        bp_1.position = pos1
    if not pos2 == bp_2.position:
        bp_2 = copy.copy(bp_2)
        bp_2.position = pos2
    
    is_dup = True if any(is_dup) else None
         
    by_genome_id_pass = defaultdict(int)
    happ_support_1 = defaultdict(list)
    happ_support_2 = defaultdict(list)
    unique_read_keys = sorted(unique_reads, key=lambda k: len(unique_reads[k]), reverse=True)
    for key, values in unique_reads.items():
        if unique_reads_pass[key]:
            by_genome_id_pass[key[0]] += len(unique_reads_pass[key])
            happ_support_1[key[0]].append(key[1])
            happ_support_2[key[0]].append(key[2])
            
    if by_genome_id_pass.values():
        is_pass = 'PASS'
        if max(by_genome_id_pass.values()) < min_reads:
            is_pass = 'FAIL'
        for keys in unique_read_keys:
            genome_id = keys[0]
            haplotype_1 = keys[1]
            haplotype_2 = keys[2]
            
            if sum(happ_support_1[genome_id]) == 3 or sum(happ_support_2[genome_id]) == 3 or sum(happ_support_1[genome_id]) == sum(happ_support_2[genome_id]) == 0:
                genotype = 'hom'
            else:
                genotype = 'het'
    
            support_reads = list(unique_reads[keys])
            supp = len(support_reads)
            
            if bp_1.ref_id == bp_2.ref_id:
                length_bp = abs(bp_1.position - bp_2.position)
                if length_bp < sv_size:
                    continue
            else:
                length_bp = 0
            prec = 1
            if not bp_1.prec or not bp_2.prec:
                prec = 0
            db_list.append(DoubleBreak(bp_1, bp_1.dir_1, bp_2, bp_2.dir_1,genome_id, haplotype_1, haplotype_2, supp, support_reads, length_bp, genotype , 'dashed'))
            db_list[-1].prec = prec
            db_list[-1].is_dup = is_dup
            db_list[-1].is_pass = is_pass
    return db_list


def match_del(double_breaks):
    DEL_THR = 20000
    CLUSTER_SIZE = 50

    clusters = defaultdict(list)
    db_list = defaultdict(list)

    for db in double_breaks:
        if db.bp_1.ref_id == db.bp_2.ref_id and db.bp_1.dir_1 == 1 and db.bp_2.dir_1 == -1 and db.length <= DEL_THR:
            clusters[db.to_string()].append(db)
            
    for cl in clusters.values():
        db = cl[0]
        db_list[db.bp_1.ref_id].append((db.bp_1.position, cl))

    clus = []
    for bp_pos in db_list.values():
        bp_pos.sort(key=lambda x:x[0])
        cluster1 = []
        for pos in bp_pos:
            if not cluster1:
                cluster1 = [pos[1]]
            elif pos[1][0].bp_1.position < cluster1[-1][0].bp_2.position < pos[1][0].bp_2.position and abs(cluster1[-1][0].length - pos[1][0].length) < CLUSTER_SIZE:
                cluster1.append(pos[1])
            else:
                if len(cluster1) > 1:
                    clus.append(cluster1)
                cluster1 = [pos[1]]
                
        if cluster1 and len(cluster1) > 1:
            clus.append(cluster1)
        
    for cl_list in clus:
        cl_0 = cl_list[0]
        cl0 = defaultdict(list)
        for db0 in cl_0:
            cl0[(db0.genome_id, db0.haplotype_1)] = db0
        for cl in cl_list[1:]:
            for db in cl:
                if (db.genome_id, db.haplotype_1) in list(cl0.keys()):
                    db0 = cl0[(db.genome_id, db.haplotype_1)]
                    db0.supp_read_ids += db.supp_read_ids
                    db0.supp = len(set(db0.supp_read_ids))
                    db.is_pass = 'FAIL_IMPREC_DEL'
                    db0.prec = 0
                else:
                    db.bp_1 = copy.copy(db.bp_1)
                    db.bp_2 = copy.copy(db.bp_2)
                    db.bp_1.position = cl_0[0].bp_1.position
                    db.bp_2.position = cl_0[0].bp_2.position
                    db.length = cl[0].length
                    cl0[(db.genome_id, db.haplotype_1)] = db
                    
    db_ls = []
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    match_breakends(double_breaks)
    for cl in clusters.values():
        if 'PASS' in [db.is_pass for db in cl]:
            db_ls += cl
    return db_ls
        
        


def match_breakends(double_breaks):
    CLUSTER_SIZE = 50
    db_list = defaultdict(list)
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)


    for cl in clusters.values():
        db = cl[0]
        db_list[db.bp_1.ref_id].append((db.bp_1.dir_1 * db.bp_1.position, 1, cl))
        db_list[db.bp_2.ref_id].append((db.bp_2.dir_1 * db.bp_2.position, 2, cl))

    clus = []
    for bp_pos in db_list.values():
        bp_pos.sort(key=lambda x:x[0])
        cluster1 = []
        for pos in bp_pos:
            if cluster1 and abs(cluster1[-1][0] - pos[0]) > CLUSTER_SIZE:
                if len(cluster1) > 1 and len(set([cl[1] for cl in cluster1])) == 2:
                    
                    clus.append(cluster1)
                cluster1 = [pos]
            else:
                cluster1.append(pos)
        if cluster1 and len(cluster1) > 1 and len(set([cl[1] for cl in cluster1])) == 2:
            clus.append(cluster1)

    for cluster1 in clus:
        conn = list()
        bp1list = [cl for cl in cluster1 if cl[1] == 1]
        for bp1 in bp1list:
            conn += bp1[2][0].bp_1.connections
        bp2list = [cl for cl in cluster1 if cl[1] == 2]
        for bp2 in bp2list:
            conn += bp2[2][0].bp_2.connections
        conn = list(set(conn))
        for cl in bp1list:
            for db in cl[2]:
                db.bp_1.connections = conn
        for cl in bp2list:
            for db in cl[2]:
                db.bp_2.connections = conn


#TODO BAM/HAPLOTYPE SPECIFIC FILTER
def double_breaks_filter(double_breaks, min_reads, control_id):

    PASS_2_FAIL_RAT = 0.5
    CONN_2_SUPP_RAT = 0.3
    CONN_2_PASS = 0.5
    CHR_CONN = 2
    COV_THR = 3
    MIN_MAPQ = 30
    NUM_HAPLOTYPES = [0,1,2]
    
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    db_list = []
    for cl in clusters.values():
        by_genome_id = defaultdict(list)
        for db in cl:
            by_genome_id[db.genome_id].append(db)
        for genome_id, dbs in by_genome_id.items():
            to_keep = []
            haplotype1 = [db.haplotype_1 for db in dbs]
            haplotype2 = [db.haplotype_2 for db in dbs]
            if 0 in haplotype1 and sum(haplotype1) > 0:
                to_keep = [i for i, hp1 in enumerate(haplotype1) if hp1]
                hp = to_keep[0]
                span_bp1 = dbs[0].bp_1.spanning_reads[genome_id][0]
                dbs[hp].bp_1.spanning_reads[genome_id][haplotype1[hp]] += span_bp1
                for db in dbs:
                    db.bp_1.spanning_reads[genome_id][0]= 0
                for db in [db for db in dbs if db.haplotype_1 == 0]:
                    db.supp = 0
                    dbs[hp].supp_read_ids += db.supp_read_ids
                dbs[hp].supp = len(set(dbs[hp].supp_read_ids))
                
            if 0 in haplotype2 and sum(haplotype2) > 0:
                if to_keep:
                    hp2 = [hp for i, hp in enumerate(haplotype2) if hp == 0 and i in to_keep]
                    if hp2:
                        hp2ls = defaultdict(list)
                        for ind in to_keep:
                            hp2ls[haplotype1[ind]].append((haplotype2[ind], ind))
                        for val in hp2ls.values():
                            hps = [c[0] for c in val]
                            if not 0 in hps:
                                continue
                            if hps == [0]:
                                new_hp = [hp for hp in haplotype2 if not hp == 0]
                                dbs[val[0][1]].haplotype_2 = new_hp[0]
                            else:
                                ind_unphased = [b for (a,b) in val if a]
                                ind_phased = [b for (a,b) in val if not a]
                                db0 = dbs[ind_unphased[0]]
                                dbs2 = dbs[ind_phased[0]]
                                dbs2.supp_read_ids += db0.supp_read_ids
                                dbs2.supp = len(set(dbs2.supp_read_ids))
                                db0.supp = 0
                                dbs2.bp_1.spanning_reads[genome_id][ind_phased[0]] += db0.bp_1.spanning_reads[genome_id][0]
                                db0.bp_1.spanning_reads[genome_id][0]= 0
                                ind0 = [i for i, db in enumerate(dbs) if db == db0]
                                to_keep.remove(ind0)
                else:
                    to_keep = [i for i, hp2 in enumerate(haplotype2) if hp2]
                    hp = to_keep[0]
                    span_bp1 = dbs[0].bp_2.spanning_reads[genome_id][0]
                    dbs[hp].bp_2.spanning_reads[genome_id][haplotype2[hp]] += span_bp1
                    for db in dbs:
                        db.bp_2.spanning_reads[genome_id][0]= 0
                    for db in [db for db in dbs if db.haplotype_2 == 0]:
                        db.supp = 0
                        dbs[hp].supp_read_ids += db.supp_read_ids
                    dbs[hp].supp = len(set(dbs[hp].supp_read_ids))
                    
            if not to_keep:
                to_keep = list(range(len(dbs)))
            db_list += [dbs[i] for i in list(set(to_keep))]
            
    for cl in clusters.values():
        db = cl[0]
        if not 'PASS' in [db.is_pass for db in cl]:
            continue
        conn_1 = db.bp_1.connections
        conn_2 = db.bp_2.connections#
        conn_pass_1 =[cn for cn in conn_1 if cn.is_pass1 == 'PASS']
        conn_pass_2 =[cn for cn in conn_2 if cn.is_pass2 == 'PASS']#
        conn_count_1 = Counter([cn.is_pass1 for cn in conn_1])
        conn_count_2 = Counter([cn.is_pass2 for cn in conn_2])#
        supp_read = sum([db.supp for db in cl])
        if conn_count_1['PASS'] < len(conn_1) * PASS_2_FAIL_RAT or conn_count_2['PASS'] < len(conn_2) * PASS_2_FAIL_RAT:
            for db1 in cl:
                db1.is_pass = 'FAIL_MAP_CONS'
            continue#
        if supp_read < conn_count_1['PASS'] * CONN_2_SUPP_RAT or supp_read < conn_count_2['PASS'] * CONN_2_SUPP_RAT:
            for db1 in cl:
                db1.is_pass = 'FAIL_CONN_CONS'
            continue#    
        conn_valid_1 = Counter([len(cn.bp_list) for cn in conn_pass_1])
        conn_valid_2 = Counter([len(cn.bp_list) for cn in conn_pass_2])
        if conn_valid_1[2] < len(conn_pass_1) * CONN_2_PASS and conn_valid_2[2] < len(conn_pass_2) * CONN_2_PASS:
            for db1 in cl:
                db1.is_pass = 'FAIL_CONN_CONS'
            continue#
        conn_ref_1 = Counter([cn.ref_id_1 for cn in conn_pass_1])
        conn_ref_2 = Counter([cn.ref_id_2 for cn in conn_pass_2])
        if len(conn_ref_1) > CHR_CONN or len(conn_ref_2) > CHR_CONN:
            for db1 in cl:
                db1.is_pass = 'FAIL_CONN_CONS'
            continue#
        qual_list = []
        for db in cl:
            qual_list += [db.bp_1.qual, db.bp_2.qual]
        vcf_qual = np.median(qual_list)
        if vcf_qual < MIN_MAPQ:
            for db in cl:
                db.is_pass = 'FAIL_MAP_CONS'
        for db in cl:
            db.vcf_qual = vcf_qual
        for db in cl:        
            if db.supp < min_reads:
                db.is_pass = 'FAIL_LOWSUPP'
                continue#
            
        has_ins = [cn.has_ins for cn in conn_pass_1 if cn in conn_pass_2]
        if any(has_ins):
            for db1 in cl:
                db1.has_ins = np.median(has_ins)
        
        #if not conn_count_1['PASS'] or not conn_count_2['PASS']:
        #    db.is_pass = 'FAIL_LOWSUPP'
        #    continue#
    
    ## report all the events not just pass    
    if  control_id:
        for cl in clusters.values():
            if not control_id in [db1.genome_id for db1 in cl]:
                haplotype1 = list(set(db1.haplotype_1 for db1 in cl))
                haplotype1 = haplotype1 if not haplotype1 == [0] else NUM_HAPLOTYPES
                haplotype2 = list(set(db1.haplotype_2 for db1 in cl))
                haplotype2 = haplotype1 if not haplotype2 == [0] else NUM_HAPLOTYPES
                span_bp1 = 0
                span_bp2 = 0
                for i in haplotype1:
                    span_bp1 += cl[0].bp_1.spanning_reads[control_id][i]
                for i in haplotype2:    
                    span_bp2 += cl[0].bp_2.spanning_reads[control_id][i]
                if span_bp1 <= COV_THR and span_bp2 <= COV_THR:
                    for db1 in cl:
                        db1.is_pass = 'FAIL_LOWCOV_NORMAL'
    
                                   
    return db_list

def add_inbetween_ins(double_breaks):
    t = 1
    new_double_breaks = []
    for db in double_breaks:
        if not db.has_ins:
            new_double_breaks.append(db)
            continue
        contig_id = 'contig_' + str(t)
        new_bp = Breakpoint(db.bp_1.ref_id, db.bp_1.position+1, 1, 60)
        new_bp.insertion_size = db.has_ins
        new_bp.contig_id = contig_id
        t += 1
        new_double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, new_bp, -1, db.genome_id, db.haplotype_1, db.haplotype_1,
                                             db.supp, db.supp_read_ids, db.has_ins, db.genotype, db.edgestyle))
        new_double_breaks[-1].has_ins = db.has_ins
        new_double_breaks.append(DoubleBreak(db.bp_2, db.direction_2, new_bp, -1, db.genome_id, db.haplotype_2, db.haplotype_1,
                                             db.supp, db.supp_read_ids, db.has_ins, db.genotype, db.edgestyle))
        new_double_breaks[-1].has_ins = db.has_ins
    return new_double_breaks

              
def extract_insertions(ins_list, clipped_clusters,ref_lengths, args):

    CLUST_LEN = 1000
    CV_THR = 0.25
    SV_LEN_THR = 0.25
    sv_len_diff = args.bp_cluster_size
    min_reads = args.bp_min_support
    min_ref_flank = args.min_ref_flank 
    sv_size = args.min_sv_size
    MIN_FULL_READ_SUPP = 2
    NUM_HAPLOTYPES = 3
    ins_clusters = []
    
    for seq, ins_pos in ins_list.items():
     
        clipped_to_remove = []
        clusters = []
        cur_cluster = []
        ins_pos.sort(key=lambda x:x.ref_end)
        clipped_clusters_seq = clipped_clusters[seq]
        clipped_clusters_seq.sort(key=lambda x:x.position)
        clipped_clusters_pos = [bp.position for bp in clipped_clusters_seq]
        for rc in ins_pos:
            if cur_cluster and rc.ref_end - cur_cluster[-1].ref_end > CLUST_LEN:
                cur_cluster.sort(key=lambda x:x.segment_length)
                cl_ins = []
                for cl1 in cur_cluster:
                    if cl_ins and abs(cl1.segment_length - cl_ins[-1].segment_length) > max(sv_len_diff, cl_ins[-1].segment_length * SV_LEN_THR):
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
            unique_reads_pass = defaultdict(set)
            unique_reads = defaultdict(set)
            unique_reads_read = defaultdict(set)
            prec = 1
            for x in cl:
                unique_reads[(x.genome_id, x.haplotype)].add(x)
                unique_reads_read[(x.genome_id, x.haplotype)].add(x.read_id)
                if x.is_pass == 'PASS':
                    unique_reads_pass[(x.genome_id, x.haplotype)].add(x)#
                    
            by_genome_id_pass = defaultdict(int)
            happ_support_1 = defaultdict(list)#
            
            for key, values in unique_reads.items():
                if unique_reads_pass[key]:
                    by_genome_id_pass[key[0]] += len(set([red.read_id for red in unique_reads_pass[key]]))
                    happ_support_1[key[0]].append(key[1])#
                    
            if not by_genome_id_pass.values() or max(by_genome_id_pass.values()) < MIN_FULL_READ_SUPP:
                continue
            
            pos_list = [x.ref_end_ori for reads in unique_reads_pass.values() for x in reads]
            pos_list.sort()
            s_len = [x.segment_length for reads in unique_reads_pass.values() for x in reads]
            ins_length = int(np.median(s_len))
            
            if ins_length * (len(s_len) + 1) < abs(pos_list[-1] - pos_list[0]):
                continue
            
            if max(pos_list) - min(pos_list) > 2 * sv_size:
                prec = 0
            elif np.std(s_len) / np.mean(s_len) > CV_THR:
                prec = 0
                
            position = int(np.median(pos_list))
            mapq = int(np.median([x.mapq for x in cl if x.is_pass == 'PASS']))
            
            ins_seq_loc = [i for i , x in enumerate(cl) if x.ins_seq and ':' not in x.ins_seq]
            if not ins_seq_loc:
                ins_seq_loc = [i for i , x in enumerate(cl) if x.ins_seq]
                
            ins_seq = cl[int(np.median(ins_seq_loc))].ins_seq
            if ins_length < sv_size:
                continue
            
            if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank:
                cl = add_clipped_end(position, clipped_clusters_pos, clipped_clusters_seq, by_genome_id_pass,
                                                         happ_support_1, unique_reads_pass)
                if cl:
                    clipped_to_remove.append(cl)
                if not by_genome_id_pass.values() or not max(by_genome_id_pass.values()) >= min_reads:
                    continue#
                for key, unique_reads in unique_reads.items():
                    #unique_reads = unique_reads_pass[key]
                    bp_1 = Breakpoint(seq, position, -1, mapq)
                    bp_1.read_ids = [x.read_id for x in unique_reads]
                    bp_1.connections = unique_reads
                    bp_3 = Breakpoint(seq, position, 1, mapq)
                    bp_3.is_insertion = True
                    bp_3.insertion_size = ins_length
                    supp = len(unique_reads_read[key])
                    supp_reads = [s.read_id for s in unique_reads]
                    genome_id = key[0]
                    if sum(happ_support_1[genome_id]) == NUM_HAPLOTYPES or sum(happ_support_1[genome_id]) == 0:
                        genotype = 'hom'
                    else:
                        genotype = 'het'
                    db_1 = DoubleBreak(bp_1, -1, bp_3, 1, genome_id, key[1], key[1], supp, supp_reads, ins_length, genotype, 'dashed')
                    db_1.prec = prec
                    db_1.ins_seq = ins_seq
                    ins_clusters.append(db_1)
        if clipped_to_remove:              
            for clipped in list(set(clipped_to_remove)):
                clipped_clusters[seq].remove(clipped)
            
    return ins_clusters


def insertion_filter(ins_list, min_reads, control_id):
    PASS_2_FAIL_RAT = 0.5
    COV_THR = 3
    NUM_HAPLOTYPES = [0,1,2]
    MIN_MAPQ = 30
    
    clusters = defaultdict(list) 
    for br in ins_list:
        clusters[br.to_string()].append(br)
    
    for cl in clusters.values():
        if not cl[0].is_pass == 'PASS':
            continue
        conn_1 = [cn for ins in cl for cn in ins.bp_1.connections]
        conn_count_1 = Counter([cn.is_pass for cn in conn_1])
        
        if not conn_count_1['PASS']:
            for ins1 in cl:
                ins1.is_pass = 'FAIL_MAP_CONS'
            continue
        
        if conn_count_1['PASS'] < len(conn_1) * PASS_2_FAIL_RAT:
            for ins1 in cl:
                ins1.is_pass = 'FAIL_MAP_CONS'
            continue
        
        qual_list = []
        for db in cl:
            qual_list += [db.bp_1.qual, db.bp_2.qual]
        vcf_qual = np.median(qual_list)
        
        if vcf_qual < MIN_MAPQ:
            for db in cl:
                db.is_pass = 'FAIL_MAP_CONS'
        for db in cl:
            db.vcf_qual = vcf_qual
                
        
        #if ins.supp < min_reads:
        #    for ins1 in cl:
        #        ins1.is_pass = 'FAIL_LOWSUPP'
        #    continue

    if control_id:
        for cl in clusters.values():
            if not control_id in [db1.genome_id for db1 in cl]:
                haplotype1 = list(set(db1.haplotype_1 for db1 in cl))
                haplotype1 = haplotype1 if not haplotype1 == [0] else NUM_HAPLOTYPES
                span_bp1 = 0
                for i in haplotype1:
                    span_bp1 += cl[0].bp_1.spanning_reads[control_id][i]
                if span_bp1 <= COV_THR:
                    for ins1 in cl:
                        ins1.is_pass = 'FAIL_LOWCOV_NORMAL'
    
    db_list = []
    for cl in clusters.values():
        by_genome_id = defaultdict(list)
        for db in cl:
            by_genome_id[db.genome_id].append(db)
        for genome_id, dbs in by_genome_id.items():
            to_keep = list(range(len(dbs)))
            haplotype1 = [db.haplotype_1 for db in dbs]
            if 0 in haplotype1 and sum(haplotype1) > 0:
                to_keep = [i for i, hp1 in enumerate(haplotype1) if hp1]
                hp = to_keep[0]
                supp = [db.supp for db in dbs if db.haplotype_1 == 0]
                span_bp1 = dbs[0].bp_1.spanning_reads[genome_id][0]
                dbs[hp].supp += sum(supp)
                dbs[hp].bp_1.spanning_reads[genome_id][haplotype1[hp]] += span_bp1
                for db in dbs:
                    db.bp_1.spanning_reads[genome_id][0] = 0
                for db in [db for db in dbs if db.haplotype_1 == 0]:
                    db.supp = 0
                    dbs[hp].supp_read_ids += db.supp_read_ids
            db_list += [dbs[i] for i in list(set(to_keep))]
    
    return db_list
    

def get_clipped_reads(segments_by_read):
    clipped_reads = defaultdict(list)
    for read in segments_by_read:
        for seg in read:
            if seg.is_clipped and seg.is_pass == 'PASS':
                clipped_reads[seg.ref_id].append(seg)
    return clipped_reads
    

def cluster_clipped_ends(clipped_reads, clust_len, min_ref_flank, ref_lengths):
    bp_list = defaultdict(list)
    QUAL = 60
    for seq, read in clipped_reads.items():
        clusters = []
        cur_cluster = []
        read.sort(key=lambda x:x.ref_end)
        for rc in read:
            if cur_cluster and rc.ref_end - cur_cluster[-1].ref_end > clust_len: 
                clusters.append(cur_cluster)
                cur_cluster = [rc]
            else:
                cur_cluster.append(rc)
        if cur_cluster:
            clusters.append(cur_cluster)
            
        for cl in clusters:
            position = int(np.median([x.ref_end_ori for x in cl]))
            if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank:
                bp = Breakpoint(seq, position, cl[0].strand, QUAL)
                bp.insertion_size = int(np.median([x.segment_length for x in cl]))
                bp.connections = cl
                bp_list[seq].append(bp)
                
    return bp_list


def add_clipped_end(position, clipped_clusters_pos, clipped_clusters_seq, by_genome_id_pass, happ_support_1, unique_reads_pass):
    
    ind = bisect.bisect_left(clipped_clusters_pos, position)
    cl = []
    MIN_SV_DIFF = 50
    if ind < len(clipped_clusters_pos)-1 and abs(clipped_clusters_pos[ind] - position) < MIN_SV_DIFF:
        cl = clipped_clusters_seq[ind]
    elif ind > 0 and abs(clipped_clusters_pos[ind - 1] - position) < MIN_SV_DIFF:
        cl = clipped_clusters_seq[ind - 1]
        
    if cl:
        cl.pos2.append(position)
        for x in cl.connections:
            if x.is_pass == 'PASS':
                unique_reads_pass[(x.genome_id,x.haplotype)].add(x)
                
        for key, values in unique_reads_pass.items():
            if unique_reads_pass[key]:
                by_genome_id_pass[key[0]] = len(unique_reads_pass[key])
                happ_support_1[key[0]].append(key[1])
    return cl

def add_breakends(double_breaks, clipped_clusters, min_reads):
    db_list = defaultdict(list)
    MIN_SV_DIFF = 100
    NUM_HAPLOTYPES = 3
    
    for db in double_breaks:
        db_list[db.bp_1.ref_id].append(db.bp_1.position)
        db_list[db.bp_2.ref_id].append(db.bp_2.position)
    t = 1
    for seq , clipped_clusters_seq in clipped_clusters.items():
        db_pos = list(set(db_list[seq]))
        for bp_1 in clipped_clusters_seq:
            position = bp_1.position
            ind = bisect.bisect_left(db_pos, position)
            cl = []
            unique_reads_pass = defaultdict(set)
            by_genome_id_pass = defaultdict(set)
            happ_support_1 = defaultdict(list)
            if ind < len(db_pos)-1 and abs(db_pos[ind] - position) > MIN_SV_DIFF:
                cl = bp_1.connections
            elif ind > 0 and abs(db_pos[ind - 1] - position) > MIN_SV_DIFF:
                cl = bp_1.connections
            
            if cl:
                for x in cl:
                    if x.is_pass == 'PASS':
                        unique_reads_pass[(x.genome_id,x.haplotype)].add(x)
                        
                for key, values in unique_reads_pass.items():
                    if unique_reads_pass[key]:
                        by_genome_id_pass[key[0]] = len(unique_reads_pass[key])
                        happ_support_1[key[0]].append(key[1])
                        
                if by_genome_id_pass.values() and max(by_genome_id_pass.values()) >= min_reads:
                    loose_end_id = 'loose_end_' + str(t)
                    t += 1
                    for key, unique_reads in unique_reads_pass.items():
                        new_bp = Breakpoint(bp_1.ref_id, bp_1.position, 1, 60) #### direction!
                        new_bp.loose_end_id = loose_end_id
                        new_bp.insertion_size = bp_1.insertion_size
                        bp_1.insertion_size = None
                        supp = len(unique_reads)
                        supp_reads = unique_reads
                        genome_id = key[0]
                        if sum(happ_support_1[genome_id]) == NUM_HAPLOTYPES:
                            genotype = 'hom'
                        else:
                            genotype = 'het'
                        double_breaks.append(DoubleBreak(bp_1, -1, new_bp, 1, genome_id, key[1], key[1], supp, supp_reads, bp_1.insertion_size, genotype, 'dashed'))
                        double_breaks[-1].sv_type = 'loose_end'

def ins_to_tra (ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size):
    
    INS_WIN = 2000 
    NUM_HAPLOTYPE = 3
    
    ins_1 = ins_list_pos[bp1.ref_id]
    strt = bisect.bisect_left(ins_1, bp1.position - INS_WIN)
    end = bisect.bisect_left(ins_1, bp1.position + INS_WIN)
    flag = False
    
    med_seg_len = int( np.quantile([min(cn.seg_len) for cn in bp2.connections],0.90))
    
    if strt == end:
        return []
    
    ins_db = ins_list[bp1.ref_id]
    clusters2 = defaultdict(list)
    for ins in ins_db[strt:end]:
        clusters2[ins.to_string()].append(ins)
    
    for ins_cl in clusters2.values():
        gen_id_1 = defaultdict(list)
        hp_list = defaultdict(list)
        ins = ins_cl[0]
        if (ins.length + min_sv_size) < abs(ins.bp_1.position - bp1.position) or (ins.length + min_sv_size) < med_seg_len:
            continue
        
        flag = True
        
        if dbs[0].bp_1.ref_id == dbs[0].bp_2.ref_id:
            svtype = 'Templated_ins'
        else:
            svtype = 'Templated_ins'
            
        for ins in ins_cl:
            hp_list[ins.genome_id].append(ins.haplotype_1)
        
        for db in dbs:
            gen_id_1[(db.genome_id, db.haplotype_1)].append(db)
            hp_list[db.genome_id].append(db.haplotype_2)
            db.sv_type = svtype
            
        for ins in ins_cl:
            ins.is_pass = 'FAIL_LONG'
            genotype = 'hom' if sum(set(hp_list[ins.genome_id])) == NUM_HAPLOTYPE else 'het'
            if gen_id_1[(ins.genome_id, ins.haplotype_1)]:
                db = gen_id_1[(ins.genome_id, ins.haplotype_1)][0]
                n_sup = list(set(ins.supp_read_ids) - set(db.supp_read_ids))
                db.genotype = genotype
                db.supp_read_ids += n_sup
                db.supp = len(set(db.supp_read_ids))
            else:
                double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, db.bp_2, db.direction_2 ,ins.genome_id, ins.haplotype_1, ins.haplotype_1, ins.supp, ins.supp_read_ids, ins.length, genotype , 'dashed'))
                double_breaks[-1].sv_type = svtype
            
    return flag

def tra_to_ins(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size):
    
    INS_WIN = 2000 
    NUM_HAPLOTYPE = 3
    total_supp_thr = 2
    
    ins_1 = ins_list_pos[bp1.ref_id]
    strt = bisect.bisect_left(ins_1, bp1.position - INS_WIN)
    end = bisect.bisect_left(ins_1, bp1.position + INS_WIN)
    flag = False
    
    med_seg_len = int( np.quantile([min(cn.seg_len) for cn in bp2.connections],0.90))
    
    if strt == end:
        return []
    
    ins_db = ins_list[bp1.ref_id]
    clusters = defaultdict(list)
    for ins in ins_db[strt:end]:
        clusters[ins.to_string()].append(ins)
    
    for ins_cl in clusters.values():
        gen_id_1 = defaultdict(list)
        hp_list = defaultdict(list)
        ins = ins_cl[0]
        total_supp = 0
        if (ins.length + min_sv_size) < abs(ins.bp_1.position - bp1.position) or (ins.length + min_sv_size) < med_seg_len:
            continue
        
        flag = True
        svtype = f"INS:{bp2.ref_id}:{bp2.position}"
            
        for ins in ins_cl:
            gen_id_1[(ins.genome_id, ins.haplotype_1)].append(ins)
            hp_list[ins.genome_id].append(ins.haplotype_1)
            ins.sv_type = svtype
        
        for db in dbs:
            hp_list[db.genome_id].append(db.haplotype_1)
        
        for db in dbs:
            genotype = 'hom' if sum(set(hp_list[db.genome_id])) == NUM_HAPLOTYPE else 'het'
            if gen_id_1[(db.genome_id, db.haplotype_1)]:
                ins = gen_id_1[(db.genome_id, db.haplotype_1)][0]
                n_sup = list(set(db.supp_read_ids) - set(ins.supp_read_ids))
                ins.supp += len(n_sup)
                ins.supp_read_ids += n_sup
                total_supp += len(n_sup)
                ins.genotype = genotype
            else:
                ins_clusters.append(DoubleBreak(ins.bp_1, ins.direction_1, ins.bp_2, ins.direction_2, db.genome_id, db.haplotype_1,  db.haplotype_1, db.supp, db.supp_read_ids, ins.length, genotype , 'dashed'))
                ins_clusters[-1].sv_type = svtype
                
        if total_supp > total_supp_thr:
            for db in dbs:
                db.is_pass = 'FAIL_LONG'
            
            
    return flag                   

def conv_tra_ins(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size, tra_vs_ins):
    if tra_vs_ins:
        return tra_to_ins(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size)
    else:
        return ins_to_tra(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size)
        
        
        
def dup_to_ins(ins_list_pos, ins_list, dbs, min_sv_size, ins_clusters, double_breaks):
    
    NUM_HAPLOTYPE = 3
    db = dbs[0]
    INS_LEN_THR = 1.2
    
    ins_1 = ins_list_pos[db.bp_1.ref_id]
    strt, end = bisect.bisect_left(ins_1, db.bp_1.position), bisect.bisect_left(ins_1, db.bp_2.position)
    
    if end - strt < 1:
        return []
    
    ins_db = ins_list[db.bp_1.ref_id]
    
    clusters = defaultdict(list)
    for ins in ins_db[strt:end]:
        clusters[ins.to_string()].append(ins)
        
    for ins_cl in clusters.values():
        gen_id_1 = defaultdict(list)
        hp_list = defaultdict(list)
        
        if db.length > ins_cl[0].length * INS_LEN_THR:
            continue
        
        for ins in ins_cl:
            hp_list[ins.genome_id].append(ins.haplotype_1)
        
        for db in dbs:
            gen_id_1[(db.genome_id, db.haplotype_1)].append(db)
            hp_list[db.genome_id].append(db.haplotype_1)
            
        for ins in ins_cl:
            ins.is_pass = 'FAIL_LONG'
            genotype = 'hom' if sum(set(hp_list[ins.genome_id])) == NUM_HAPLOTYPE else 'het'
            if gen_id_1[(ins.genome_id, ins.haplotype_1)]:
                db = gen_id_1[(ins.genome_id, ins.haplotype_1)][0]
                n_sup = list(set(ins.supp_read_ids) - set(db.supp_read_ids))
                db.supp += len(n_sup)
                db.supp_read_ids += n_sup
                db.genotype = genotype
            else:
                double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, db.bp_2, db.direction_2 ,ins.genome_id, ins.haplotype_1, ins.haplotype_1, ins.supp, ins.supp_read_ids, ins.length, genotype , 'dashed'))
                double_breaks[-1].is_dup = True
  
def match_long_ins(ins_clusters, double_breaks, min_sv_size, tra_vs_ins):

    DEL_THR = 100000
    ins_list = defaultdict(list)
    ins_list_pos = defaultdict(list)
    for ins in ins_clusters:
        ins_list_pos[ins.bp_1.ref_id].append(ins.bp_1.position)
        ins_list[ins.bp_1.ref_id].append(ins)
        
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
        
    for dbs in clusters.values():
        db = dbs[0]
        if db.bp_1.ref_id == db.bp_2.ref_id and  db.direction_1 == db.direction_2:
            continue
        if db.bp_1.ref_id == db.bp_2.ref_id and db.direction_1 > 0 and db.direction_2 < 0 and db.bp_2.position - db.bp_1.position < DEL_THR:
            continue
        if db.is_dup:
            dup_to_ins(ins_list_pos, ins_list, dbs, min_sv_size, ins_clusters, double_breaks)
        else:
            ins_to_remove_tra = conv_tra_ins(ins_list_pos, ins_list, db.bp_1, db.bp_2, db.direction_1, dbs, ins_clusters, double_breaks, min_sv_size, tra_vs_ins)
            if not ins_to_remove_tra:
                conv_tra_ins(ins_list_pos, ins_list, db.bp_2, db.bp_1, db.direction_2, dbs, ins_clusters, double_breaks, min_sv_size, tra_vs_ins)
                
def add_insseq(double_breaks, segments_by_read, bam_files, thread_pool):
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
        
    db_ls = [cl for cl in clusters.values() if cl[0].has_ins]
    tasks = [(segments_by_read, bam_files, cl) for cl in db_ls]
    parsing_results = None
    parsing_results = thread_pool.starmap(get_insseq, tasks)
    
    for [cl,ins_seq] in parsing_results:
        for db in cl:
            db.ins_seq = ins_seq
        
        
def get_insseq (segments_by_read, bam_files, cl):
    db = cl[0]
    cn = [conn for conn in db.bp_1.connections if conn in db.bp_2.connections]
    conn = cn[np.argmin(np.array([c.has_ins for c in cn]) - db.has_ins)]
    read_id = conn.read_id
    for seg in segments_by_read:
        if seg and seg[0].read_id == read_id:
            break
    for s in seg:
        if s.is_primary:
            break
    bam_file = bam_files[s.genome_id]
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(s.ref_id, s.ref_start-5, s.ref_end+1,  multiple_iterators=True):
        if aln.query_name == read_id:
            break
    ins_seq = aln.query_sequence[conn.ins_pos:(conn.ins_pos + int(db.has_ins))]
    return(cl, ins_seq)

def calc_vaf(db_list):
    for db1 in db_list.values():
        db = db1[0]
        span_bp1 = sum(db.bp_1.spanning_reads[db.genome_id])
        span_bp2 = sum(db.bp_2.spanning_reads[db.genome_id])
        DV = 0
        DR = int(np.mean([span_bp1, span_bp2])) if not db.bp_2.is_insertion else span_bp1
        for db in db1:
            DR1 = int(np.median([db.bp_1.spanning_reads[(db.genome_id, db.haplotype_1)], db.bp_2.spanning_reads[(db.genome_id, db.haplotype_2)]])) if not db.bp_2.is_insertion else db.bp_1.spanning_reads[(db.genome_id, db.haplotype_1)]
            DV1 = db.supp
            db.hvaf =  DV1 / (DV1 + DR1) if DV1 > 0 else 0
            DV += DV1
        vaf =  DV / (DV + DR) if DV > 0 else 0
        for db in db1:
            db.DR = DR
            db.DV = DV
            db.vaf = vaf

def add_mut_type(db_list, control_id, control_vaf):
    mut_type = 'germline'#
    sample_ids = list(db_list.keys())
    if not sample_ids == [control_id] and (not control_id in sample_ids or db_list[control_id][0].vaf < control_vaf):
        mut_type = 'somatic'
    for genome_id, db1 in db_list.items():
        if genome_id == control_id:
            for db in db1:
                db.mut_type = 'germline'
        else:
            pass_list = [db.is_pass for db in db1]
            if 'FAIL_LOWCOV_OTHER' in pass_list and not 'PASS' in pass_list:
                mut_type = 'germline'
            for db in db1:
                db.mut_type = mut_type
        
    
def annotate_mut_type(double_breaks, control_id, control_vaf, vaf_thr):
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    for db_clust in clusters.values():
        vaf_pass = 'FAIL'
        db_list = defaultdict(list)
        
        for db in db_clust:
            db_list[db.genome_id].append(db)
        
        calc_vaf(db_list)
        vaf_list = [db.vaf for db in db_clust]
        if max(vaf_list) >= vaf_thr:
            vaf_pass = 'PASS'
            
        for db in db_clust:
            db.vaf_pass = vaf_pass
        
        if control_id:
            add_mut_type(db_list, control_id, control_vaf)


def add_sv_type(double_breaks):
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    t = 0
    for db_clust in clusters.values():
        sv_type = get_sv_type(db_clust[0])
        sv_id = 'severus_' + sv_type + str(t)
        t +=1
        
        for db in db_clust:    
            db.vcf_sv_type = sv_type
            db.vcf_id = sv_id
    
def get_sv_type(db):
    if db.bp_2.loose_end_id:
        db.sv_type = 'loose_end'
        return 'BND'
    if db.sv_type == 'tandem_duplication':
        return 'DUP'
    if db.bp_2.is_insertion or db.bp_1.is_insertion:
        return 'INS'
    if not db.bp_1.ref_id == db.bp_2.ref_id:
        return 'BND'
    if db.direction_1 == db.direction_2:
        return 'INV'
    if db.bp_1.dir_1 > 0:
        return 'DEL'
    else:
        return 'BND'        

def filter_germline_db(double_breaks):
    db_list = defaultdict(list) 
    for db in double_breaks:
        db_list['germline'].append(db)
        if db.mut_type == 'somatic':
            db_list['somatic'].append(db)
    return db_list

def filter_fail_double_db(double_breaks, output_all, vaf_thr, min_sv_size, coverage_histograms, segments_by_read, bam_files, thread_pool, ins_seq):
    db_list = []
    if not output_all:
        for db in double_breaks:
            if db.is_pass == 'PASS' and db.vaf_pass == 'PASS':
                db_list.append(db)
    else:
        db_list = double_breaks
    cluster_db(db_list, coverage_histograms, min_sv_size)
    add_sv_type(double_breaks)
    if ins_seq:
        add_insseq(db_list, segments_by_read, bam_files, thread_pool)
    db_list = filter_germline_db(db_list)
    return db_list

def get_phasingblocks(hb_vcf):
    MIN_BLOCK_LEN = 10000
    MIN_SNP = 10
    
    vcf = pysam.VariantFile(hb_vcf)
    haplotype_blocks = defaultdict(list)
    startpoint_list = defaultdict(list)
    endpoint_list = defaultdict(list)
    switch_points = defaultdict(list)
    id_list = defaultdict(list)

    for var in vcf:
        if 'PS' in var.samples.items()[0][1].items()[-1]:
            haplotype_blocks[(var.chrom, var.samples.items()[0][1]['PS'])].append(var.pos)

    phased_lengths = []
    for (chr_id, block_name), coords in haplotype_blocks.items():
        if max(coords) - min(coords) > MIN_BLOCK_LEN and len(coords) >= MIN_SNP:
            startpoint_list[chr_id].append(min(coords))
            endpoint_list[chr_id].append(max(coords))
            phased_lengths.append(max(coords) - min(coords))
            id_list[chr_id].append(block_name)

    for chr_id, end_list in endpoint_list.items():
        start_list = startpoint_list[chr_id]
        switch_points[chr_id] = [(a + b) // 2 for a, b in zip(start_list[1:], end_list[:-1])]

    total_phased = sum(phased_lengths)
    _l50, n50 = _calc_nx(phased_lengths, total_phased, 0.50) 
    logger.info(f"\tTotal phased length: {total_phased}")
    logger.info(f"\tPhase blocks N50: {n50}")

    return (switch_points, id_list)


def add_phaseset_id(double_breaks, hb_points):
    (switch_points, id_list) = hb_points
    for db in double_breaks:
        ind_1 = bisect.bisect_left(switch_points[db.bp_1.ref_id], db.bp_1.position)
        ind_2 = bisect.bisect_left(switch_points[db.bp_2.ref_id], db.bp_2.position)
        id1 = id_list[db.bp_1.ref_id][ind_1] if ind_1 < len(id_list[db.bp_1.ref_id]) else 0
        id2 = id_list[db.bp_2.ref_id][ind_2] if ind_2 < len(id_list[db.bp_2.ref_id]) else 0
        db.phaseset_id = (id1, id2)
    
def segment_coverage(histograms, genome_id, ref_id, ref_start, ref_end, haplotype):
    hist_start = ref_start // COV_WINDOW
    hist_end = ref_end // COV_WINDOW
    cov_list = histograms[(genome_id, haplotype, ref_id)][hist_start : hist_end + 1]
    if not cov_list:
        return 0
    return int(np.median(cov_list))

def get_segments_coverage(db_segments, coverage_histograms, max_genomic_length):
    genomic_segments = defaultdict(list)
    phasing_segments = defaultdict(list)
    
    for db, segments in db_segments.items():
        for (genome_id, seg_ref, seg_start, seg_end, seg_hp, sw_point) in segments:
            
            if seg_end - seg_start > max_genomic_length:
                if seg_start in [db.bp_1.position,  db.bp_2.position]:
                    seg_end = seg_start + max_genomic_length
                else: 
                    seg_start = seg_end - max_genomic_length
            if seg_ref == 'INS':
                coverage = db.supp
                gs = GenomicSegment(genome_id, seg_hp, db.bp_1.ref_id, db.bp_1.position, db.bp_1.position,
                                                       coverage, seg_end - seg_start)
                gs.is_insertion = True
            else:
                coverage = segment_coverage(coverage_histograms, genome_id, seg_ref, seg_start, seg_end, seg_hp)
            
                if db.is_dup and seg_end - seg_start == 0:
                    continue
                    
                gs = GenomicSegment(genome_id, seg_hp, seg_ref, seg_start, seg_end,
                                                       coverage, seg_end - seg_start)
            genomic_segments[db].append(gs)
            
            if sw_point:
                phasing_segments[sw_point].append(gs)
                
    return (genomic_segments, phasing_segments)

def get_ref_adj(genomic_segments, ref_adj):
    by_genome = defaultdict(list)
    
    for db, gs_ls in genomic_segments.items():
        for gs in gs_ls:
            if gs.pos1 in ref_adj[gs.ref_id]:
                by_genome[(gs.genome_id, gs.haplotype, gs.pos1)].append(gs)
            elif gs.pos2 in ref_adj[gs.ref_id]:
                by_genome[(gs.genome_id, gs.haplotype, gs.pos2)].append(gs)
    
    adj_segments = []
    for key, gs_ls in by_genome.items():
        pos = key[2]
        pos1 = []
        pos2 = []
        for gs in gs_ls:
            if gs.pos1 == pos:
                pos1.append(gs)
            else:
                pos2.append(gs)
        for ps in pos1:
            for ps2 in pos2:
                adj_segments.append((ps, ps2))
                
    return adj_segments
 
def get_genomic_segments(double_breaks, coverage_histograms, hb_vcf, key_type, ref_lengths, min_ref_flank, max_genomic_length):
    switch_points = defaultdict(list)
    hb_points = []
    if hb_vcf:
        hb_points = get_phasingblocks(hb_vcf)
        switch_points = hb_points[0]
    
    bp1_list = defaultdict(list)
    bp2_list = defaultdict(list)
    single_bp = defaultdict(list)
    db_segments = defaultdict(list)
    ref_adj = defaultdict(list)
    
    for double_bp in double_breaks:
        bp1_list[(double_bp.genome_id, double_bp.haplotype_1, double_bp.bp_1.ref_id)].append(double_bp)
        
        if not double_bp.bp_2.is_insertion:
            bp2_list[(double_bp.genome_id, double_bp.haplotype_2, double_bp.bp_2.ref_id)].append(double_bp)
            
        single_bp[(double_bp.genome_id, double_bp.haplotype_1, double_bp.bp_1.ref_id)].append((double_bp.bp_1.position, double_bp.direction_1, ''))  
        single_bp[(double_bp.genome_id, double_bp.haplotype_2, double_bp.bp_2.ref_id)].append((double_bp.bp_2.position, double_bp.direction_2, ''))
        
        if double_bp.is_dup:
            single_bp[(double_bp.genome_id, double_bp.haplotype_1, double_bp.bp_1.ref_id)].append((double_bp.bp_1.position, 1, 'dup'))
            single_bp[(double_bp.genome_id, double_bp.haplotype_2, double_bp.bp_2.ref_id)].append((double_bp.bp_2.position, -1, 'dup'))

    for (genome_name, haplotype_name, ref_name), s_bp in single_bp.items():
        for ps in s_bp:
            if ps[2]:
                ref_adj[(ref_name)].append(ps[0])
        sw_p = [(sw,0) for sw in switch_points[ref_name]]
        sw_p += [(min_ref_flank, 0), (ref_lengths[ref_name] - min_ref_flank,0 )]
        s_bp += sw_p
        s_bp = list(set(s_bp))
        s_bp.sort(key=lambda x: (x[0],x[1]))
        pos_ls = [s[0] for s in s_bp if not s[1] == -1]
        neg_ls = [s[0] for s in s_bp if not s[1] == 1]
        sw_p = [s[0] for s in s_bp if s[1] == 0]
        bp1 = bp1_list[(genome_name, haplotype_name, ref_name)]
        for db in bp1:
            sw_point = ''
            if db.bp_2.is_insertion:
                db_segments[db].append((genome_name, 'INS', 0, db.length, haplotype_name,sw_point))
                ind = bisect.bisect_left(pos_ls,db.bp_1.position)
                if pos_ls[ind -1] in sw_p:
                    sw_point = (ref_name, pos_ls[ind -1])
                db_segments[db].append((genome_name, ref_name, pos_ls[ind -1], db.bp_1.position,  haplotype_name, sw_point))
                sw_point = ''
                ind2 = bisect.bisect_left(neg_ls,db.bp_1.position)
                if neg_ls[ind2 +1] in sw_p:
                    sw_point = (ref_name, neg_ls[ind2 +1])
                db_segments[db].append((genome_name, ref_name, db.bp_1.position, neg_ls[ind2 +1], haplotype_name, sw_point))
            elif db.bp_1.dir_1 == -1:
                ind = bisect.bisect_right(pos_ls,db.bp_1.position)
                if not db.is_dup and db.bp_1.position == pos_ls[ind]:
                    ind = ind +1
                if pos_ls[ind] in sw_p:
                    sw_point = (ref_name, pos_ls[ind])
                db_segments[db].append((genome_name, ref_name, db.bp_1.position, pos_ls[ind], haplotype_name, sw_point))
            elif db.bp_1.dir_1 == 1:
                ind = bisect.bisect_left(neg_ls,db.bp_1.position)
                if neg_ls[ind -1] == db.bp_1.position:
                    ind = ind -1
                if neg_ls[ind - 1] in sw_p:
                    sw_point = (ref_name, neg_ls[ind -1])
                db_segments[db].append((genome_name, ref_name, neg_ls[ind -1], db.bp_1.position, haplotype_name, sw_point))
        bp2 = bp2_list[(genome_name, haplotype_name, ref_name)]
        for db in bp2:
            sw_point = ''
            if db.bp_2.dir_1 == -1:
                ind = bisect.bisect_right(pos_ls,db.bp_2.position)
                if db.bp_2.position == pos_ls[ind]:
                    ind = ind +1
                if pos_ls[ind] in sw_p:
                    sw_point = (ref_name, pos_ls[ind])
                db_segments[db].append((genome_name, ref_name, db.bp_2.position, pos_ls[ind], haplotype_name, sw_point))
            elif db.bp_2.dir_1 == 1:
                ind = bisect.bisect_left(neg_ls,db.bp_2.position)
                if not db.is_dup and neg_ls[ind -1] == db.bp_2.position:
                    ind = ind -1
                if neg_ls[ind - 1] in sw_p:
                    sw_point = (ref_name, neg_ls[ind -1])
                db_segments[db].append((genome_name, ref_name, neg_ls[ind -1], db.bp_2.position, haplotype_name, sw_point))
                
    (genomic_segments, phasing_segments) = get_segments_coverage(db_segments, coverage_histograms, max_genomic_length)
    
    adj_segments = get_ref_adj(genomic_segments, ref_adj)
    
    if key_type == 'germline' and hb_points:
        add_phaseset_id(double_breaks, hb_points)
        
    return (genomic_segments, phasing_segments, adj_segments)

def get_insertionreads(segments_by_read):
    ins_list_all = defaultdict(list)
    for read in segments_by_read:
        for seg in read:
            if seg.is_insertion:
                ins_list_all[seg.ref_id].append(seg)
    return ins_list_all    

def get_splitreads(segments_by_read):
    split_reads = []
    for read in segments_by_read:
        split = [seg for seg in read if not seg.is_insertion and not seg.is_clipped]
        if len(split)>1:
            split_reads.append(split)
    return split_reads

def resolve_overlaps(segments_by_read, min_ovlp_len):
    """
    Some supplementary alignments may be overlapping (e.g. in case of inversions with flanking repeat).
    This function checks if the overlap has ok structe, trims and outputs non-overlapping alignments
    """
    
    def _get_ovlp(seglist_1, seglist_2):
        
        alg_strt = [seglist_1[-1].align_start, seglist_2[-1].align_start]
        alg_end = [seglist_1[-1].read_end, seglist_2[-1].read_end]
        
        if alg_end[0] - alg_strt[1] > min_ovlp_len:
            return  alg_end[0] - alg_strt[1] + 1
        else:
            return 0
        
    def _get_full_ovlp(seglist_1, seglist_2):
        
        alg_strt = [seglist_1[-1].align_start, seglist_2[-1].align_start]
        alg_end = [seglist_1[-1].read_end, seglist_2[-1].read_end]
        
        max_ovlp_len = [alg_end[0] - alg_strt[0], alg_end[1] - alg_strt[1]]
        if alg_end[0] - alg_strt[1] == min(max_ovlp_len):
            seg = [seglist_1, seglist_2]
            return seg[max_ovlp_len.index(min(max_ovlp_len))]
        else:
            return False
            
    def _update_ovlp_seg(seglist_1, seglist_2, left_ovlp):
        seg_to_remove = []
        
        seg2 = seglist_1 if seglist_2[0].strand == 1 else seglist_2
        
        for seg in seg2:
            seg.align_start += left_ovlp
            if seg.read_start >= seg.align_start:
                continue
            if seg.read_end < seg.align_start:
                seg_to_remove.append(seg)
            else:
                seg.read_start += left_ovlp
                seg.segment_length = seg.read_end - seg.read_start
                seg.ref_start = seg.ref_start + left_ovlp
                seg.ref_start_ori = seg.ref_start
        
        return seg_to_remove

    for read in segments_by_read:
        read_segments = [seg for seg in read if not seg.is_insertion and not seg.is_clipped]
        if len(read_segments) < 2:
            continue
    
        segs_to_remove =[]
        cur_cluster = []
        clusters = []
        read_segments.sort(key=lambda s:(s.align_start, s.read_start))
        
        for seg in read_segments:
            if cur_cluster and not seg.align_start == cur_cluster[-1].align_start: 
                clusters.append(cur_cluster)
                cur_cluster = [seg]
            else:
                cur_cluster.append(seg)
        if cur_cluster:
            clusters.append(cur_cluster)
            
        if len(clusters) < 2:
            continue
        
        for i in range(len(clusters)):
            left_ovlp = 0
            seg = False
            
            if i > 0 and clusters[i-1][0].ref_id == clusters[i][0].ref_id:
                seg_to_remove = _get_full_ovlp(clusters[i-1], clusters[i])
                if seg_to_remove:
                    segs_to_remove += seg_to_remove
                    continue
                left_ovlp = _get_ovlp(clusters[i-1], clusters[i])
                
            if left_ovlp > 0:
                seg_to_remove = _update_ovlp_seg(clusters[i-1], clusters[i], left_ovlp)
                segs_to_remove += seg_to_remove
                
        for seg in list(set(segs_to_remove)):
            read.remove(seg)
               
def cluster_db(double_breaks2, coverage_histograms, min_sv_size):
    clusters = defaultdict(list)
    for br in double_breaks2:
        br.subgraph_id = []
        if br.bp_1.is_insertion:
            continue
        clusters[br.to_string()].append(br)
    
    conn_duplications(clusters, coverage_histograms)
    conn_inter(clusters)
         
def cluster_inversions(double_breaks, coverage_histograms, min_sv_size):
    by_genome = defaultdict(list)
    for db in double_breaks:
        by_genome[db.genome_id].append(db)
    
    for db_list in by_genome.values():
        complex_inv(db_list, coverage_histograms, min_sv_size)

def complex_inv(double_breaks, coverage_histograms, min_sv_size):
    inv_list = defaultdict(list)
    
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    foldback_inv(clusters, coverage_histograms)
    
    for cl in clusters.values():
        db = cl[0]
        if db.bp_1.dir_1 == db.bp_2.dir_1:
            inv_list[db.bp_1.ref_id].append(cl)
        if not db.bp_1.ref_id == db.bp_2.ref_id:
            inv_list[db.bp_2.ref_id].append(cl)
    
    for seq, inv in inv_list.items():
        pos_list = [(cl[0].bp_1.position, cl[0].bp_1.dir_1, cl[0].haplotype_1, cl[0].phaseset_id[0], cl) for cl in inv if cl[0].bp_1.ref_id == seq and not cl[0].sv_type] + \
            [(cl[0].bp_2.position, cl[0].bp_1.dir_1,cl[0].haplotype_2, cl[0].phaseset_id[1], cl) for cl in inv if cl[0].bp_2.ref_id == seq and not cl[0].sv_type]
        
        pos_list.sort(key=lambda s:(s[0], -s[1]))
        
        neg_lis = [ind for ind, (_,_,_,_,db) in enumerate(pos_list) if db[0].bp_1.dir_1 == -1]
        
        pairs = []
        pairs_pos = []
        for i in neg_lis:
            if i > len(pos_list)-2:
                continue
            black_list = []
            if len(pos_list[i][4]) == 1 and not pos_list[i][2] == 0:
                black_list = [ind for ind, (_,dir1, hp, phase_id, _) in enumerate(pos_list) if ind > i and phase_id == pos_list[i][3] and not hp == pos_list[i][2]]
            for j in range (i + 1, len(pos_list)):
                if not j in black_list:
                    break
            if pos_list[j][4][0].bp_1.dir_1 == 1:
                pairs.append((pos_list[i][4],pos_list[j][4]))
                pairs_pos.append((pos_list[i][0],pos_list[j][0]))
        
        node_list = defaultdict(int)
        G = nx.Graph()
        t = 0
        for (db1, db2) in pairs:
            if not db1[0] in node_list.keys():
                node_list[db1[0]] = t
                G.add_node(t, _cluster = db1)
                t += 1
            rn = node_list[db1[0]]
            if not db2[0] in node_list.keys():
                node_list[db2[0]] = t
                G.add_node(t, _cluster = db2)
                t += 1
            ln = node_list[db2[0]]
            G.add_edge(rn,ln)
        
        
        for i, c in enumerate(nx.connected_components(G)):
            cl_list = []
            for n_id in c:   
                cl_list.append(G.nodes[n_id]['_cluster'])
                    
            sv_type = 'complex_inv'
            if len(cl_list) == 2:
                db1 = cl_list[0][0] if cl_list[0][0].bp_1.dir_1 == 1 else cl_list[1][0]
                db2 = cl_list[1][0] if cl_list[0][0].bp_1.dir_1 == 1 else cl_list[0][0]
                
                if not db1.bp_1.ref_id == db1.bp_2.ref_id or not db2.bp_1.ref_id == db2.bp_2.ref_id:
                    sv_type = 'inv_tra'
                    
                else:
                    comb = [(db1.bp_1.position,db1.bp_1.dir_1), (db1.bp_2.position,db1.bp_2.dir_1), (db2.bp_1.position,db2.bp_1.dir_1), (db2.bp_2.position,db2.bp_2.dir_1)]
                    comb.sort(key=lambda s:s[0])
                    dir_list = [c[1] for c in comb]
                    if db1.bp_1.position - min_sv_size < db2.bp_1.position < db1.bp_2.position < db2.bp_2.position + min_sv_size:
                        if (db1.bp_2.position - db2.bp_1.position) / (db2.bp_2.position - db1.bp_1.position) > 0.90:
                            sv_type = 'reciprocal_inv'
                        else:
                            sv_type = 'reciprocal_inv_del'
                    elif db2.bp_1.position - min_sv_size < db1.bp_1.position < db2.bp_2.position + min_sv_size:
                        sv_type = 'dup_inv_segment'
                    
                    elif dir_list == [-1, -1, 1, 1] or dir_list == [-1, 1, 1, -1]:
                        sv_type = 'dup_inv_segment'
                        
                    elif dir_list == [1, -1, -1, 1]:
                        sv_type = 'inv_inserion'
            for cl in cl_list:
                for db in cl:
                    db.cluster_id = 'INV' + str(i)
                    db.sv_type = sv_type
            
                
def foldback_inv(clusters, coverage_histograms):
    inv_list_pos = defaultdict(list)
    inv_list_neg = defaultdict(list)
    FOLD_BACK_THR = 0.5
    FOLD_BACK_DIST_THR = 50000
    t = 0
    
    for cl in clusters.values():
        db = cl[0]
        if not db.bp_2.ref_id == db.bp_1.ref_id:
            continue
        if db.direction_1 == db.direction_2 == 1:
            inv_list_pos[db.bp_1.ref_id].append(cl)
            
        elif db.direction_1 == db.direction_2 == -1:
            inv_list_neg[db.bp_1.ref_id].append(cl)
            
    for seq, pos_lis in inv_list_pos.items():
        neg_lis = inv_list_neg[seq]
        foldback_pairs = defaultdict(list)
        for ind, cl in enumerate(pos_lis):
            db = cl[0]
            if db.bp_2.position - db.bp_1.position > FOLD_BACK_DIST_THR:
                continue
            pos1 = db.bp_1.position // 500
            pos2 = db.bp_2.position // 500
            max_pos = len(coverage_histograms[(db.genome_id, db.haplotype_2, db.bp_2.ref_id)])
            max_pos1 = len(coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)])
            cov1 = coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos1-1:min(pos1+2,max_pos1)]
            cov2 = coverage_histograms[(db.genome_id, db.haplotype_2, db.bp_2.ref_id)][min(pos2+2,max_pos)]
            thr = int(db.supp * FOLD_BACK_THR)
            if cov1[0] > cov1[2] + thr and  cov1[2] > cov2 + thr:
                foldback_pairs['pos'].append(ind)
                for db in pos_lis[ind]:
                    db.sv_type = 'foldback'
    
        for ind, cl in enumerate(neg_lis):
            db = cl[0]
            if db.bp_2.position - db.bp_1.position > FOLD_BACK_DIST_THR:
                continue
            pos1 = db.bp_1.position // 500
            pos2 = db.bp_2.position // 500
            max_pos = len(coverage_histograms[(db.genome_id, db.haplotype_2, db.bp_2.ref_id)])
            max_pos1 = len(coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)])
            cov1 = coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos1-1:min(pos1+2,max_pos1)]
            cov2 = coverage_histograms[(db.genome_id, db.haplotype_2, db.bp_2.ref_id)][min(pos2+2,max_pos)]
            thr = int(db.supp * FOLD_BACK_THR)
            if cov2 > cov1[2] + thr and  cov1[2] > cov1[0] + thr:
                foldback_pairs['neg'].append(ind)
                for db in neg_lis[ind]:
                    db.sv_type = 'foldback'
        
        if len(foldback_pairs) == 2:
            t += 1
            for ind in foldback_pairs['pos']:
                for db in pos_lis[ind]:
                    db.sv_type = 'BFB_foldback'
                    db.cluster_id = 'BFB' + str(t)
                    
            for ind in foldback_pairs['neg']:
                for db in neg_lis[ind]:
                    db.sv_type = 'BFB_foldback'
                    db.cluster_id = 'BFB' + str(t)
            
        
            
def conn_duplications(clusters, coverage_histograms):        
    
    DUP_COV_THR = 0.5
    COV_WINDOW = 500
    for ind, cl in enumerate(clusters.values()):
        db = cl[0]
        is_dup = False
        if not db.bp_1.ref_id == db.bp_2.ref_id or db.bp_2.is_insertion:
            continue
        
        if db.is_dup:
            db.sv_type = 'tandem_duplication'
            continue
        
        if db.direction_1 == -1 and db.direction_2 == 1:
            pos1 = db.bp_1.position // COV_WINDOW
            pos2 = db.bp_2.position // COV_WINDOW
            
            max_pos = len(coverage_histograms[(db.genome_id, db.haplotype_2, db.bp_2.ref_id)])
            cov1 = coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][max(pos1-5,0):pos1]
            cov1 += coverage_histograms[(db.genome_id, db.haplotype_2, db.bp_2.ref_id)][pos2:min(pos2+5,max_pos)]
            cov1 = int(np.median(cov1))
            cov3 = int(np.median(coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos1:pos2+1]))
            
            if cov3 > cov1 + db.supp * DUP_COV_THR:
                is_dup = True
            svtype = 'tandem_duplication' if is_dup else 'Templated_ins'
            for db in cl:
                db.is_dup = is_dup
                db.sv_type = svtype
            
def conn_inter(clusters):
    intra_chr_cl = defaultdict(list)
    DEL_THR = 1000000
    t = 0
    
    for cl in clusters.values():
        db = cl[0]
        if db.bp_1.ref_id == db.bp_2.ref_id:
            continue
        intra_chr_cl[(db.bp_1.ref_id, db.bp_2.ref_id)].append(cl)
        
    for key, intra_chr in intra_chr_cl.items():
        if len(intra_chr) < 2:
            continue
        
        ind_list = []
        for ind, cl in enumerate(intra_chr):
            if ind in ind_list:
                continue
            db = cl[0]
            dir1, dir2 = -1 * db.direction_1, -1 * db.direction_2
            for ind2, cl2 in enumerate(intra_chr[ind+1:]):
                db2 = cl2[0]
                if db2.direction_1 == dir1 and db2.direction_2 == dir2:
                    dist = max([abs(db.bp_1.position - db2.bp_1.position), abs(db.bp_2.position - db2.bp_2.position)]) 
                    if dist > DEL_THR:
                        continue
                    ind_list.append(ind2)
                    t += 1
                    for db in cl:
                        db.sv_type = 'Templated_ins'
                        db.cluster_id = 'TIC' + str(t)
                        
                    for db2 in cl2:
                        db2.sv_type = 'Templated_ins'
                        db.cluster_id = 'TIC' + str(t)

def write_alignments(allsegments, outpath):
    aln_dump_stream = open(outpath, "w")
    for read in allsegments:
        for seg in read:
            if seg.is_insertion or seg.is_clipped:
                continue
            aln_dump_stream.write(str(seg) + "\n")
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
                def_array.append(('',0,0,0))
                loc[(tag,k)] = t
                header += '_'.join([tag, str(k),'PASS/FAIL,'])
                header += '_'.join([tag, str(k),'support,'])
                header += '_'.join([tag, str(k),'spanning_1,'])
                header += '_'.join([tag, str(k),'spanning_2,'])
                t += 1
    else:
        for tag in genome_tags:
            def_array.append(('',0,0,0))
            loc[(tag,0)] = t
            header += '_'.join([tag, 'PASS/FAIL,'])
            header += '_'.join([tag, 'support,'])
            header += '_'.join([tag, 'spanning_1,'])
            header += '_'.join([tag, 'spanning_2,'])
            t += 1
            
    summary_csv = defaultdict(list)
    for br in double_breaks:
        if not summary_csv[br.to_string()]:
            summary_csv[br.to_string()] = def_array[:]
            idd=(br.genome_id, br.haplotype_1)
            summary_csv[br.to_string()][loc[idd]] = (br.is_pass, br.supp, br.bp_1.spanning_reads[br.genome_id][br.haplotype_1], br.bp_2.spanning_reads[br.genome_id][br.haplotype_2])            
        else:
            idd=(br.genome_id, br.haplotype_1)
            summary_csv[br.to_string()][loc[idd]] = (br.is_pass, br.supp, br.bp_1.spanning_reads[br.genome_id][br.haplotype_1], br.bp_2.spanning_reads[br.genome_id][br.haplotype_2])
            
    out_stream.write(header + "\n")
    for key,values in summary_csv.items():
        bp_array=[]
        for k in values:
            bp_array.append(str(k[0]))
            bp_array.append(str(k[1]))
            bp_array.append(str(k[2]))
            bp_array.append(str(k[3]))
        bp_to_write = ','.join([key, ','.join(bp_array)])
        out_stream.write(bp_to_write)
        out_stream.write("\n")

def output_readids(double_breaks, genome_ids, out_stream):
    header = ','.join(['#SV_ID'] + genome_ids)
    out_stream.write(header + "\n")
    clusters = defaultdict(list)
    for db in double_breaks:
        clusters[db.vcf_id].append(db)
    
    for key, cl in clusters.items():
        by_genome_id = defaultdict(list)
        by_genome_id_db = defaultdict(list)
        for db in cl:
            by_genome_id[db.genome_id] += db.supp_read_ids
            by_genome_id_db[db.genome_id].append(db)
            line = [key]
        for genome_id in genome_ids:
            if by_genome_id[genome_id]:
                read_ids = list(set(by_genome_id[genome_id]))
                db = by_genome_id_db[genome_id][0]
                line.append('\"' + '; '.join(read_ids) + '\"')
            else:
                line.append('\"\"')
        line = ','.join(line)
        out_stream.write(line)
        out_stream.write("\n")
                            
def call_breakpoints(segments_by_read, ref_lengths, coverage_histograms, bam_files, genome_ids, control_id, thread_pool, args):
    
    if args.write_alignments:
        outpath_alignments = os.path.join(args.out_dir, "read_alignments")
        write_alignments(segments_by_read, outpath_alignments)
        
    
    if args.resolve_overlaps:
        logger.info('Resolving overlaps')
        #resolve_overlaps(segments_by_read,  args.sv_size)
        
    logger.info('Extracting split alignments')
    split_reads = get_splitreads(segments_by_read)
    ins_list_all = get_insertionreads(segments_by_read)
    cont_id  = list(control_id)[0] if control_id else '' 
    
    logger.info('Extracting clipped reads')

    clipped_clusters = []
    extract_clipped_end(segments_by_read)
    clipped_reads = get_clipped_reads(segments_by_read)
    clipped_clusters = cluster_clipped_ends(clipped_reads, args.bp_cluster_size,args.min_ref_flank, ref_lengths)
    
    logger.info('Clustering unmapped insertions')
    ins_clusters = extract_insertions(ins_list_all, clipped_clusters, ref_lengths, args)
    
    logger.info('Starting breakpoint detection')
    double_breaks = get_breakpoints(split_reads, ref_lengths, args)
    
    logger.info('Starting match_long_ins')
    match_long_ins(ins_clusters, double_breaks, args.min_sv_size, args.tra_to_ins)
    
    logger.info('Starting compute_bp_coverage')
    get_coverage_parallel(bam_files, genome_ids, thread_pool, args.min_mapping_quality, double_breaks)
    get_coverage_parallel(bam_files, genome_ids, thread_pool, args.min_mapping_quality, ins_clusters)
    
    logger.info('Filtering breakpoints')
    double_breaks = double_breaks_filter(double_breaks, args.bp_min_support, cont_id)
    double_breaks.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position, b.direction_1))
    
    ins_clusters = insertion_filter(ins_clusters, args.bp_min_support, cont_id)
    ins_clusters.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position))
   
    double_breaks +=  ins_clusters
    
    if args.inbetween_ins:
        add_breakends(double_breaks, clipped_clusters, args.bp_min_support)
    
    logger.info('annotate_mut_type')
    
    annotate_mut_type(double_breaks, cont_id, args.control_vaf, args.vaf_thr)
        
    logger.info('Writing breakpoints')
    output_breaks(double_breaks, genome_ids, args.phase_vcf, open(os.path.join(args.out_dir,"breakpoints_double.csv"), "w"))
    
    double_breaks = filter_fail_double_db(double_breaks, args.output_all, args.vaf_thr, args.min_sv_size, coverage_histograms, segments_by_read, bam_files, thread_pool, args.ins_seq)
    
    if args.output_read_ids:
            output_readids(double_breaks['germline'], genome_ids, open(os.path.join(args.out_dir,"read_ids.csv"), "w"))
    
    return double_breaks 
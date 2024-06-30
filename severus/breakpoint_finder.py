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
from severus.resolve_vntr import read_vntr_file

logger = logging.getLogger()


MAX_LOWMAPQ_READS = 10
MIN_SEGMENT_LENGTH = 100
MIN_SEGMENT_OVERLAP = 100
MAX_SEGMENT_OVERLAP = 500
MAX_CONNECTION= 1000
MAX_UNALIGNED_LEN = 500
COV_WINDOW = 1000
CHUNK_SIZE = 10000000
        
class Breakpoint(object):
    __slots__ = ("ref_id", "position","dir_1", "spanning_reads", "connections", 'prec',
                 "read_ids", "pos2", 'id', "is_insertion", "insertion_size", "qual")
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
        self.prec = 1

    def fancy_name(self):
        if not self.is_insertion:
            return self.unique_name()
        return f"INS:{self.insertion_size}"

    def unique_name(self):
        if self.is_insertion:
            return f"INS:{self.ref_id}:{self.position}"
        else:
            return f"{self.ref_id}:{self.position}"

    def coord_tuple(self):
        sign = '-' if self.dir_1 == -1 else '+'
        return (self.ref_id, self.position, sign)


class DoubleBreak(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "genome_id","haplotype_1",'haplotype_2',"supp",'supp_read_ids',
                 'length','genotype','edgestyle', 'is_pass', 'ins_seq', 'mut_type', 'is_dup', 'has_ins','subgraph_id', 'sv_type','tra_pos',
                 'DR', 'DV', 'hvaf', 'vaf', 'prec', 'phaseset_id', 'cluster_id', 'gr_id', 'vcf_id', 'vcf_sv_type','vaf_pass', 'vcf_qual', 'haplotypes', 'is_single', 'vntr')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, genome_id, haplotype_1, haplotype_2, 
                 supp, supp_read_ids, length):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
        self.direction_2 = direction_2
        self.genome_id = genome_id
        self.haplotype_1 = haplotype_1
        self.haplotype_2 = haplotype_2
        self.haplotypes= []
        self.supp = supp
        self.supp_read_ids = supp_read_ids
        self.length = length
        self.genotype = ''
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
        self.gr_id = 0
        self.cluster_id = 0
        self.vcf_id = None
        self.vcf_sv_type = None
        self.vaf_pass = None
        self.vcf_qual = None
        self.is_single = None
        self.vntr = None
        self.tra_pos = None
        
    def to_string(self):
        strand_1 = "+" if self.direction_1 > 0 else "-"
        strand_2 = "+" if self.direction_2 > 0 else "-"
        label_1 = "{0}{1}:{2}".format(strand_1, self.bp_1.ref_id, self.bp_1.position)
        if self.bp_2.is_insertion:
            label_1 = "{0}:{1}".format(self.bp_1.ref_id, self.bp_1.position)
            label_2 = "{0}:{1}".format('INS', self.length)
        elif self.is_single:
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
        elif self.is_single:
            label_2 = "{0}:{1}".format(self.bp_2.ref_id, self.bp_2.position)
        else:
            label_2 = "{0}{1}:{2}".format(strand_2, self.bp_2.ref_id, self.bp_2.position)
            if label_2[1:] < label_1[1:]:
                label_1, label_2 = label_2, label_1
        bp_name = f"{label_1}\t{label_2}"
        return bp_name
        

class GenomicSegment(object):
    __slots__ = "genome_id","haplotype", "ref_id", 'dir1', "pos1", 'dir2', "pos2", "coverage", "length_bp", "is_insertion", "total_coverage"
    def __init__(self, genome_id, haplotype, ref_id, pos1, pos2, coverage, total_coverage,length_bp):
        self.genome_id = genome_id
        self.haplotype = haplotype
        self.ref_id = ref_id
        self.dir1 = -1
        self.dir2 = 1
        self.pos1 = pos1
        self.pos2 = pos2
        self.coverage = coverage
        self.total_coverage = total_coverage
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

        
def get_pos(segls, bp_dir):
    dirls = ['right', 'left']
    s1 = segls[bp_dir]
    if not segls[1].read_start > segls[0].read_start:
        return s1.get_pos(dirls[bp_dir - 1])
    else:
        return s1.get_pos(dirls[bp_dir])
           
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
    single_bps = []
    
    seq_breakpoints_l = defaultdict(list)
    seq_breakpoints_r = defaultdict(list)
    
    def _add_double(s1, s2):
        ref_bp_1 = s1.get_pos("right")[0]
        ref_bp_2 = s2.get_pos("left")[0]
        if ref_bp_1 > ref_bp_2:
            s1, s2 = s2, s1
        rc = (s1,s2)
        seq_breakpoints_r[s1.ref_id].append(rc)
        seq_breakpoints_l[s2.ref_id].append(rc)
        
    for read_segments in split_reads:
        read_segments.sort(key=lambda x:(x.align_start, x.read_start))
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            if s2.read_start - s1.read_end < MAX_SEGMENT_DIST:
                _add_double(s1, s2)
                
    all_breaks = []
    for seq, bp_pos in seq_breakpoints_r.items():
        bps = cluster_bp(seq, bp_pos, clust_len, min_ref_flank, ref_lengths, min_reads,0)
        if bps:
            all_breaks += bps
            
    for seq, bp_pos in seq_breakpoints_l.items():  
        bps = cluster_bp(seq, bp_pos, clust_len, min_ref_flank, ref_lengths, min_reads,1)
        if bps:
            all_breaks += bps
    
    conn_list = defaultdict(list)
    for bp in all_breaks:
        for conn in bp.connections:
            conn_list[conn].append(bp)
       
    matched_bp, bp_ls = match_breaks(conn_list)
    if args.single_bp:
        single_bps = get_single_bps(bp_ls)
    double_breaks=[]
    for (bp_1 , bp_2), cl in matched_bp.items():
        db = get_double_breaks(bp_1, bp_2, cl, sv_size, min_reads,bp_ls)
        if db:
            double_breaks += db
    double_breaks = match_del(double_breaks)
    
    return double_breaks, single_bps

def get_single_bps(bp_ls):
    single_bps = []
    for bp, rcls in bp_ls.items():
        if 2 not in rcls:
            single_bps.append(bp)
    return single_bps
    
def cluster_bp(seq, bp_pos, clust_len, min_ref_flank, ref_lengths, min_reads, bp_dir):
    clusters = []
    cur_cluster = []
    bp_list = []
    min_supp = 2
    bp_pos.sort(key=lambda bp: (get_pos(bp, bp_dir)[2], get_pos(bp, bp_dir)[0]))
    for rc in bp_pos:
        if cur_cluster and get_pos(rc,bp_dir)[0] - get_pos(cur_cluster[-1], bp_dir)[0] > clust_len:
            if len(cur_cluster) >= min_supp:
                clusters.append(cur_cluster)
            cur_cluster = [rc]
        else:
            cur_cluster.append(rc)
    if cur_cluster and len(cur_cluster) >= min_supp:
        clusters.append(cur_cluster)
        
    for cl in clusters:
        unique_reads = set()
        read_ids = []
        connections =[]
        position_arr = []
        qual_arr = []
        
        for rc in cl:
            x = rc[bp_dir]
            unique_reads.add((x.read_id, (x.genome_id,x.haplotype)))
            read_ids.append(x.read_id)
            connections.append(rc)
            if x.is_pass == 'PASS':
                position_arr.append(get_pos(rc, bp_dir)[1])
                qual_arr.append(x.mapq)
            
        by_genome_id = defaultdict(int)
        for read in unique_reads:
            by_genome_id[read[1]] += 1
        if not position_arr:
            continue
        position = int(np.median(position_arr))
        prec = 1
        if max(position_arr) - min(position_arr) > 2 * clust_len:
            prec = 0
        qual = int(np.median(qual_arr))
        sign  = get_pos(rc, bp_dir)[2]
        if position >= min_ref_flank and position <= ref_lengths[seq] - min_ref_flank:
            bp = Breakpoint(seq, position, sign, qual)
            bp.connections = connections
            bp.read_ids = read_ids
            bp.prec = prec
            bp_list.append(bp)
                
    return bp_list
            
def match_breaks(conn_list):
    matched_bp = defaultdict(list)
    bp_ls = defaultdict(list)
    for rc, bp_list in conn_list.items():
        if len(bp_list) < 2:
            bp_ls[bp_list[0]].append(1)
        elif len(bp_list) == 2:
            bp_list.sort(key=lambda bp: bp.position)
            matched_bp[(bp_list[0], bp_list[1])].append(rc)
            bp_ls[bp_list[0]].append(2)
            bp_ls[bp_list[1]].append(2)
        
    return matched_bp , bp_ls
        
def get_double_breaks(bp_1, bp_2, cl, sv_size, min_reads, bp_ls):  
    unique_reads = defaultdict(set)
    unique_reads_pass = defaultdict(set)
    db_list = []
    is_dup = None
    
    CONN_2_PASS = 0.5
    
    conn_valid_1 = Counter(bp_ls[bp_1])
    conn_valid_2= Counter(bp_ls[bp_2])
    conn_pass_1 = len([cn for cn in cl if cn[0].is_pass == 'PASS'])
    conn_pass_2 = len([cn for cn in cl if cn[1].is_pass == 'PASS'])
    
    for x,y in cl:
        unique_reads[(x.genome_id,x.haplotype,y.haplotype)].add(x.read_id)
        if get_pos((x,y), 0)[2]== -1 and get_pos((x,y), 1)[2] == 1 and x.ref_id == y.ref_id and x.ref_start <= x.ref_start <= y.ref_end <= x.ref_end:
            is_dup = True
        if x.is_pass == 'PASS' and y.is_pass == 'PASS':
            unique_reads_pass[(x.genome_id,x.haplotype,y.haplotype)].add(x.read_id)
    
    pos1 = int(np.median([get_pos(c, 0)[1] for c in cl]))
    pos2 = int(np.median([get_pos(c,1)[1] for c in cl]))
    if not pos1 == bp_1.position:
        bp_1 = copy.copy(bp_1)
        bp_1.position = pos1
    if not pos2 == bp_2.position:
        bp_2 = copy.copy(bp_2)
        bp_2.position = pos2
         
    by_genome_id_pass = defaultdict(int)
    for key in unique_reads_pass.keys():
        by_genome_id_pass[key[0]] += len(unique_reads_pass[key])
            
    if by_genome_id_pass.values():
        is_pass = 'PASS'
        if max(by_genome_id_pass.values()) < min_reads:
            is_pass = 'FAIL'
            
        if conn_valid_1[2] < conn_pass_1 * CONN_2_PASS and conn_valid_2[2] < conn_pass_2 * CONN_2_PASS:
            is_pass = 'FAIL_CONN_CONS'
            
        for keys in unique_reads.keys():
            genome_id = keys[0]
            haplotype_1 = keys[1]
            haplotype_2 = keys[2]
            
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
            db_list.append(DoubleBreak(bp_1, bp_1.dir_1, bp_2, bp_2.dir_1,genome_id, haplotype_1, haplotype_2, supp, support_reads, length_bp))
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


def check_db(cl, conn_1, ind):

    PASS_2_FAIL_RAT = 0.5
    CONN_2_SUPP_RAT = 0.3
    CHR_CONN = 2
    MIN_MAPQ = 30
    db = cl[0]
    if db.bp_1.ref_id == db.bp_2.ref_id and db.bp_1.position > db.bp_2.position:
        return 'FAIL_MAP_CONS'
    conn_pass_1 =[cn for cn in conn_1 if cn[ind].is_pass == 'PASS']
    conn_count_1 = Counter([cn[ind].is_pass for cn in conn_1])
    supp_read = sum([db.supp for db in cl])
    
    if conn_count_1['PASS'] < len(conn_1) * PASS_2_FAIL_RAT:
        return 'FAIL_MAP_CONS'
    if supp_read < conn_count_1['PASS'] * CONN_2_SUPP_RAT:
        return 'FAIL_CONN_CONS'
    conn_ref_1 = Counter([cn[ind].ref_id for cn in conn_pass_1])
    if len(conn_ref_1) > CHR_CONN:
        return 'FAIL_CONN_CONS'
    qual_list = []
    for db in cl:
        qual_list += [db.bp_1.qual, db.bp_2.qual]
    vcf_qual = np.median(qual_list)
    for db in cl:
        db.vcf_qual = vcf_qual
    if vcf_qual < MIN_MAPQ:
        return 'FAIL_MAP_CONS'
    
    else:
        return ''
    
        
def double_breaks_filter(double_breaks, single_bps, min_reads, control_id, resolve_overlaps, sv_size):

    COV_THR = 3
    NUM_HAPLOTYPES = [0,1,2]
    MAX_STD = 25
    MIN_DIST = 50000
    
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
            
    for cl in clusters.values():
        db = cl[0]
        if not 'PASS' in [db.is_pass for db in cl]:
            continue
        conn_1 = db.bp_1.connections
        conn_2 = db.bp_2.connections#
        fail1 = check_db(cl, conn_1, 0)
        fail2 = check_db(cl, conn_2, 1)
        if fail1 or fail2 and 0 < db.length < MIN_DIST:
            if fail1:
                for db1 in cl:
                    db1.is_pass = fail1
            else:
                for db1 in cl:
                    db1.is_pass = fail1
            
        elif not fail1 and fail2:
            db.bp_1.pos2 = db.bp_2.unique_name()
            for db1 in cl:
                db1.bp_2 = db.bp_1
                db1.is_single = True
                single_bps.append(db1)
        elif not fail2 and fail1:
            db.bp_2.pos2 = db.bp_1.unique_name()
            for db1 in cl:
                db1.bp_1 = db.bp_2
                db1.is_single = True
                single_bps.append(db1)
        elif fail1 and fail2:
            for db1 in cl:
                db1.is_pass = fail1
                
        if not fail1 or not fail2:
            conn_ins = [cn for cn in conn_1 if cn in conn_2 and cn[0].is_pass == 'PASS' and cn[1].is_pass == 'PASS']
            has_ins = []
            for c in conn_ins:
                s1,s2 = sorted(c, key=lambda x:(x.align_start, x.read_start))
                dist = s2.read_start - s1.read_end
                has_ins.append(dist)
            has_ins_len = int(np.median(has_ins))
            if has_ins_len >= sv_size and np.std(has_ins) < MAX_STD:
                cn = conn_ins[np.argmin(np.array(has_ins) - db.has_ins)]
                s1,s2 = sorted(cn, key=lambda x:(x.align_start, x.read_start))
                for db1 in cl:
                    db1.has_ins = (s1, s2, has_ins_len)
    
    if resolve_overlaps:
        resolve_ovlp(clusters) 
       
    if  control_id:
        for cl in clusters.values():
            if cl[0].is_single:
                continue
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
    
    double_breaks = [db for db in double_breaks if not db.is_single]
    
    match_haplotypes(double_breaks)
    return double_breaks


def resolve_ovlp(clusters):
    MIN_SIZE = 50
    for cl in clusters.values():
        db = cl[0]
        if not db.is_pass == 'PASS' or db.is_single:
            continue
        conn_1 = [c for c in db.bp_1.connections if c in db.bp_2.connections]
        overlap = []
        for cn in conn_1:
            s1,s2 = sorted(cn, key=lambda x:(x.align_start, x.read_start))
            ovlp = 0
            if s1.read_end > s2.read_start:
                ovlp = s1.read_end - s2.read_start
                overlap.append(ovlp)
        if not overlap:
            continue
        ovlp = int(np.median(overlap))
        if ovlp < MIN_SIZE:
            continue
        mm1 = np.median([c[0].mismatch_rate for c in conn_1])
        mm2 = np.median([c[1].mismatch_rate for c in conn_1])
        if mm1 > mm2:
            ovlp = ovlp if db.direction_1 == -1 else -ovlp
            db.bp_1.position += ovlp
        else:
            ovlp = ovlp if db.direction_2 == -1 else -ovlp
            db.bp_2.position += ovlp
        db.length = abs(db.bp_2.position - db.bp_1.position) if db.length else 0

              
def extract_insertions(ins_list, clipped_clusters,ref_lengths, args):

    CLUST_LEN = 1000
    CV_THR = 0.25
    SV_LEN_THR = 0.25
    sv_len_diff = args.bp_cluster_size
    min_reads = args.bp_min_support
    min_ref_flank = args.min_ref_flank 
    sv_size = args.min_sv_size
    MIN_FULL_READ_SUPP = 2
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
            
            for key, values in unique_reads.items():
                if unique_reads_pass[key]:
                    by_genome_id_pass[key[0]] += len(set([red.read_id for red in unique_reads_pass[key]]))
                    
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
                cl = add_clipped_end(ins_length, position, clipped_clusters_pos, clipped_clusters_seq, by_genome_id_pass,unique_reads_pass)
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
                    
                    db_1 = DoubleBreak(bp_1, -1, bp_3, 1, genome_id, key[1], key[1], supp, supp_reads, ins_length)
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
    
    match_haplotypes(ins_list)

    

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

def get_single_bp(single_bp, clipped_clusters, double_breaks, bp_min_support, cont_id,min_ref_flank, ref_lengths):
    MIN_DIST = 100
    PASS_2_FAIL = 0.4
    filtered_sbp = []
    
    for bps in clipped_clusters.values():
        for bp in bps:
            bp.connections = [(c,c) for c in bp.connections]
            
    for sbp in single_bp:
        clipped_clusters[sbp.ref_id].append(sbp)
    
    db_pos = defaultdict(list)
    for db in double_breaks:
        db_pos[db.bp_1.ref_id].append(db.bp_1.position)
        db_pos[db.bp_2.ref_id].append(db.bp_2.position)
    
    for seq, sbp_ls in clipped_clusters.items():
        posls = db_pos[seq]
        posls += [min_ref_flank, ref_lengths[seq] - min_ref_flank]
        posls = list(set(posls))
        posls.sort()
        for s_bp in sbp_ls:
            ind = bisect.bisect_left(posls, s_bp.position)
            ind = ind if not ind == len(posls) else ind-1
            ind = ind if not ind == 0 else 1
            if abs(posls[ind] - s_bp.position) < MIN_DIST or abs(posls[ind-1] - s_bp.position) < MIN_DIST:
                continue
            conn = [c[0].is_pass for c in s_bp.connections] + [c[1].is_pass for c in s_bp.connections if not c[0] == c[1]]
            conn_count = Counter(conn)
            if conn_count['PASS'] < max([len(conn) * PASS_2_FAIL, bp_min_support]):
                continue
            cl = s_bp.connections
            unique_reads = defaultdict(set)
            unique_reads_pass = defaultdict(set)
            if len(cl) < bp_min_support:
                continue
            for x in cl:
                hp = x[0].haplotype if x[0].is_pass == 'PASS' else x[1].haplotype
                unique_reads[(x[0].genome_id,hp)].add(x[0].read_id)
                if x[0].is_pass == 'PASS' or x[1].is_pass == 'PASS':
                    unique_reads_pass[(x[0].genome_id,hp)].add(x[0].read_id)
            by_genome_id_pass = defaultdict(int)
            unique_read_keys = sorted(unique_reads, key=lambda k: len(unique_reads[k]), reverse=True)
            for key, values in unique_reads.items():
                if unique_reads_pass[key]:
                    by_genome_id_pass[key[0]] += len(unique_reads_pass[key])
            if by_genome_id_pass.values():
                if max(by_genome_id_pass.values()) < bp_min_support:
                    continue
                for keys in unique_read_keys:
                    genome_id = keys[0]
                    haplotype_1 = keys[1]
                    
                    support_reads = list(unique_reads[keys])
                    supp = len(support_reads)
                    length_bp = 0
                    prec = 1
                    if not s_bp.prec:
                        prec = 0
                    filtered_sbp.append(DoubleBreak(s_bp, s_bp.dir_1, s_bp, s_bp.dir_1,genome_id, haplotype_1, haplotype_1, supp, support_reads, length_bp))
                    filtered_sbp[-1].prec = prec
                    filtered_sbp[-1].is_pass = 'PASS'
                    filtered_sbp[-1].vcf_sv_type = 'BND'
                    filtered_sbp[-1].is_single = True
                    filtered_sbp[-1].vcf_qual = s_bp.qual
                    
    return filtered_sbp

def filter_single_bp(single_bps, cont_id, control_vaf, vaf_thr):
    sbp_list = []
    QUAL_THR = 30
    VAF_THR = 0.25
    match_haplotypes(single_bps)
    annotate_mut_type(single_bps, cont_id, control_vaf, VAF_THR)
    for sbp in single_bps:
        if sbp.vaf_pass == 'PASS' and sbp.vcf_qual > QUAL_THR and sbp.is_pass == 'PASS':
            sbp_list.append(sbp)
                        
    clusters = defaultdict(list) 
    for br in sbp_list:
        clusters[br.to_string()].append(br)
    sv_id = 'severus_sBND'
    t = 0
    for cl in clusters.values():
        for db in cl:
            db.vcf_id = sv_id + str(t)
            db.vcf_sv_type = 'BND'
            db.length = 0
            t+=1
    return sbp_list

def add_clipped_end(ins_length,position, clipped_clusters_pos, clipped_clusters_seq, by_genome_id_pass, unique_reads_pass):
    
    ind = bisect.bisect_left(clipped_clusters_pos, position)
    cl = []
    sv_diff = min(250, ins_length)
    if ind < len(clipped_clusters_pos)-1 and abs(clipped_clusters_pos[ind] - position) < sv_diff:
        cl = clipped_clusters_seq[ind]
    elif ind > 0 and abs(clipped_clusters_pos[ind - 1] - position) < sv_diff:
        cl = clipped_clusters_seq[ind - 1]
        
    if cl:
        cl.pos2.append(position)
        for x in cl.connections:
            if x.is_pass == 'PASS':
                unique_reads_pass[(x.genome_id,x.haplotype)].add(x)
                
        for key, values in unique_reads_pass.items():
            if unique_reads_pass[key]:
                by_genome_id_pass[key[0]] = len(unique_reads_pass[key])
                
    return cl

def ins_to_tra (ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size,slen):
    
    INS_WIN = 2000
    MIN_DIFF = 50
    
    ins_1 = ins_list_pos[bp1.ref_id]
    strt = bisect.bisect_left(ins_1, bp1.position - INS_WIN)
    end = bisect.bisect_left(ins_1, bp1.position + INS_WIN)
    flag = False
    
    MAX_DEPTH = 1000
    if len(bp1.connections) > MAX_DEPTH or len(bp2.connections) > MAX_DEPTH:
        med_seg_len = int( np.quantile([cn[slen].segment_length for cn in bp1.connections[:MAX_DEPTH]],0.90))
    else:
        med_seg_len = int( np.quantile([cn[slen].segment_length for cn in bp2.connections if cn in bp1.connections],0.90))
    
    #seglen = 0
    #med_seg_len = int(np.quantile([cn.seg_len[slen] for cn in bp2.connections if cn in bp1.connections],0.90))
    #seg_len = sorted([cn.seg_len[slen] for cn in bp2.connections if cn in bp1.connections])
    
    #if seg_len[-1] - seg_len[0] < MIN_DIFF:
    #    seglen = med_seg_len
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
        #if seglen and abs(ins.length - seglen) > MIN_DIFF:
        #    continue
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
            if gen_id_1[(ins.genome_id, ins.haplotype_1)]:
                db = gen_id_1[(ins.genome_id, ins.haplotype_1)][0]
                n_sup = list(set(ins.supp_read_ids) - set(db.supp_read_ids))
                db.supp_read_ids += n_sup
                db.supp = len(set(db.supp_read_ids))
            else:
                double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, db.bp_2, db.direction_2 ,ins.genome_id, ins.haplotype_1, ins.haplotype_1, ins.supp, ins.supp_read_ids, ins.length))
                double_breaks[-1].sv_type = svtype
            
    return flag

def tra_to_ins(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size, slen):
   
    INS_WIN = 2000
    total_supp_thr = 2
    MAX_DEPTH = 10000
    #MIN_DIFF = 50
    
    ins_1 = ins_list_pos[bp1.ref_id]
    strt = bisect.bisect_left(ins_1, bp1.position - INS_WIN)
    end = bisect.bisect_left(ins_1, bp1.position + INS_WIN)
    flag = False
    #seglen = 0
    db = dbs[0]
    if len(bp1.connections) > MAX_DEPTH or len(bp2.connections) > MAX_DEPTH:
        med_seg_len = int( np.quantile([cn[slen].segment_length for cn in bp1.connections[:MAX_DEPTH]],0.90))
    else:
        med_seg_len = int( np.quantile([cn[slen].segment_length for cn in bp2.connections if cn in bp1.connections],0.90))
    
    #seg_len = sorted([cn.seg_len[slen] for cn in bp2.connections if cn in bp1.connections])
    #if seg_len[-1] - seg_len[0] < MIN_DIFF:
    #    seglen = med_seg_len
    
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
        #if seglen and abs(ins.length - seglen) > MIN_DIFF:
        #    continue
        flag = True
        if bp2.dir_1 == -1:
            tra_pos = bp2.ref_id + ':' + str(bp2.position)+ '-'  + str(bp2.position + ins.length)
        else:
            tra_pos = bp2.ref_id + ':'+ str(bp2.position - ins.length) +  '-' + str(bp2.position)
            
        for ins in ins_cl:
            gen_id_1[(ins.genome_id, ins.haplotype_1)].append(ins)
            hp_list[ins.genome_id].append(ins.haplotype_1)
            ins.tra_pos = tra_pos
        
        for db in dbs:
            hp_list[db.genome_id].append(db.haplotype_1)
        
        for db in dbs:
            if gen_id_1[(db.genome_id, db.haplotype_1)]:
                ins = gen_id_1[(db.genome_id, db.haplotype_1)][0]
                n_sup = list(set(db.supp_read_ids) - set(ins.supp_read_ids))
                ins.supp += len(n_sup)
                ins.supp_read_ids += n_sup
                total_supp += len(n_sup)
            else:
                ins_clusters.append(DoubleBreak(ins.bp_1, ins.direction_1, ins.bp_2, ins.direction_2, db.genome_id, db.haplotype_1,  db.haplotype_1, db.supp, db.supp_read_ids, ins.length))
                ins_clusters[-1].tra_pos = tra_pos
                
        if total_supp > total_supp_thr:
            for db in dbs:
                db.is_pass = 'FAIL_LONG'
            
    return flag                   

def conv_tra_ins(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size, tra_vs_ins, slen):
    if tra_vs_ins:
        return tra_to_ins(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size, slen)
    else:
        return ins_to_tra(ins_list_pos, ins_list, bp1, bp2, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size, slen)
        
        
def dup_to_ins(ins_list_pos, ins_list, dbs, min_sv_size, ins_clusters, double_breaks):
    
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
            if gen_id_1[(ins.genome_id, ins.haplotype_1)]:
                db = gen_id_1[(ins.genome_id, ins.haplotype_1)][0]
                n_sup = list(set(ins.supp_read_ids) - set(db.supp_read_ids))
                db.supp += len(n_sup)
                db.supp_read_ids += n_sup
            else:
                double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, db.bp_2, db.direction_2 ,ins.genome_id, ins.haplotype_1, ins.haplotype_1, ins.supp, ins.supp_read_ids, ins.length))
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
            ins_to_remove_tra = conv_tra_ins(ins_list_pos, ins_list, db.bp_1, db.bp_2, db.direction_1, dbs, ins_clusters, double_breaks, min_sv_size, tra_vs_ins, 1)
            if not ins_to_remove_tra:
                conv_tra_ins(ins_list_pos, ins_list, db.bp_2, db.bp_1, db.direction_2, dbs, ins_clusters, double_breaks, min_sv_size, tra_vs_ins,0)
                
def add_insseq(double_breaks2, segments_by_read, bam_files, thread_pool):

    clusters = defaultdict(list) 
    for br in double_breaks2:
        clusters[br.to_string()].append(br)
    dbls = defaultdict(list)
    for cl in clusters.values():
        if not cl[0].has_ins:
            continue
        db = cl[0]
        s1,s2,ins_len = db.has_ins
        dbls[s1.read_id] = [cl, (s1.read_end, s2.read_start), ins_len]
    pos_ls = defaultdict(list)
    for seg in segments_by_read:
        if seg and seg[0].read_id in dbls.keys():
            for s in seg:
                if s.is_primary:
                    pos_ls[(s.ref_id,s.ref_start//CHUNK_SIZE, s.genome_id)].append((s.read_id, dbls[s.read_id][2], dbls[s.read_id][1]))
                    break
    tasks = [(bam_files[key[2]], key[0], key[1], val) for key, val in pos_ls.items()]
    parsing_results = None
    parsing_results = thread_pool.starmap(get_insseq, tasks)
    for res in parsing_results:
        for read_id, ins_seq in res:
            cl = dbls[read_id][0]
            for db in cl:
                db.ins_seq = ins_seq
                
def get_insseq(bam_file,ref_id, pos, val):
    ins_seq = []
    read_ids = [v[0] for v in val]
    aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(ref_id, pos * CHUNK_SIZE, (pos+1) * CHUNK_SIZE,  multiple_iterators=True):
        if aln.query_name in read_ids and not aln.is_supplementary and not aln.is_secondary and not aln.is_unmapped:
            st_pos,end_pos = [v[2] for v in val if v[0] == aln.query_name][0]
            if aln.is_reverse:
                st_pos, end_pos = aln.query_length - end_pos , aln.query_length - st_pos
            ins_seq.append((aln.query_name, aln.query_sequence[st_pos:end_pos]))
    return ins_seq

def calc_vaf(db_list):
    for db1 in db_list.values():
        db = db1[0]
        span_bp1 = sum(db.bp_1.spanning_reads[db.genome_id])
        span_bp2 = sum(db.bp_2.spanning_reads[db.genome_id])
        DV = 0
        DR = int(np.mean([span_bp1, span_bp2])) if not db.bp_2.is_insertion else span_bp1
        for db in db1:
            if db.is_pass == 'FAIL_MERGED_HP':
                continue
            DR1 = int(np.median([db.bp_1.spanning_reads[db.genome_id][db.haplotype_1], db.bp_2.spanning_reads[db.genome_id][db.haplotype_2]])) if not db.bp_2.is_insertion else db.bp_1.spanning_reads[db.genome_id][db.haplotype_1]
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
    if not sample_ids == [control_id] and (not control_id in sample_ids or db_list[control_id][0].vaf <= control_vaf):
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
        

def calc_gentype(db_list):
    for dbb in db_list.values():
        hp1 = [db.haplotype_1 for db in dbb if db.is_pass == 'PASS']
        hp2 = [db.haplotype_2 for db in dbb if db.is_pass == 'PASS']
        if hp1 and hp2:
            gentype1 = 'hom' if (sum(hp1) == 3 or sum(hp1) == 0) and (sum(hp2) == 3 or sum(hp2) == 0) else 'het'
            
            for db in dbb:
                db.genotype = gentype1
        
    
def annotate_mut_type(double_breaks, control_id, control_vaf, vaf_thr):
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    for db_clust in clusters.values():
        vaf_pass = 'FAIL'
        db_list = defaultdict(list)
        
        for db in db_clust:
            db_list[db.genome_id].append(db)
        
        calc_gentype(db_list)
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
    if db.is_single:
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

def filter_fail_double_db(double_breaks, single_bps, coverage_histograms, segments_by_read, bam_files, thread_pool, args):

    min_sv_size = args.min_sv_size
    ins_seq = args.ins_seq
    single_bp = args.single_bp
    
    db_list = []
    for db in double_breaks:
        if db.is_pass == 'PASS' and db.vaf_pass == 'PASS' and db.vcf_qual:
            db_list.append(db)
        
    cluster_db(db_list, coverage_histograms, min_sv_size)
    add_sv_type(db_list)
    
    if ins_seq:
        add_insseq(db_list, segments_by_read, bam_files, thread_pool)
        
    if single_bp:
        db_list += single_bps
    
    db_list = filter_germline_db(db_list)
    return db_list

def check_vntr(db, vntr_list):
    tr_reg = vntr_list[db.bp_1.ref_id]
    if tr_reg:
        strt = bisect.bisect_right(tr_reg[2], db.bp_1.position)
        end = bisect.bisect_left(tr_reg[3], db.bp_2.position)
        if strt - end == 1:
            return strt
        
def add_vntr_annot(double_breaks, args):
    vntr_list = read_vntr_file(args.vntr_file)
    add_sv_type(double_breaks)
    clusters = defaultdict(list)
    vntr_clusters = defaultdict(list)
    MERGE_THR = 150
    MERGE_RAT = 0.25
    for br in double_breaks:
        if br.vcf_sv_type == 'DEL' or br.vcf_sv_type == 'INS':
            clusters[br.to_string()].append(br)
    
    for cl in clusters.values():
        vntr = check_vntr(cl[0], vntr_list)
        if vntr:
            vntr_clusters[(cl[0].bp_1.ref_id, vntr)].append(cl)
            for db in cl:
                db.vntr = True
            if cl[0].vcf_sv_type == 'DEL':
                pass_conn = [1 for (a,b) in cl[0].bp_1.connections if not a.is_pass == b.is_pass == 'vntr_only']
                if len(pass_conn) < 2:
                    for db in cl:
                        db.is_pass = 'FAIL_VNTR'
            else:
                pass_conn = [1 for a in cl[0].bp_1.connections if not a.is_pass == 'vntr_only']
                if len(pass_conn) < 2:
                    for db in cl:
                        db.is_pass = 'FAIL_VNTR'
    
    for vntr_cls in vntr_clusters.values():
        if len(vntr_cls) > 3:
            for cl in vntr_cls:
                for db in cl:
                    db.is_pass = 'FAIL_COMPLEX_VNTR'
        elif len(vntr_cls) > 1:
            dbs_type = defaultdict(list)
            for cl in vntr_cls:
                dbs_type[cl[0].vcf_sv_type].append(cl)
            for key, del_dbs in dbs_type.items():
                if len(del_dbs) > 1:
                    del_dbs.sort(key=lambda s: s[0].length)
                    clusters = []
                    cur_cluster = []
                    for cl in del_dbs:
                        if cur_cluster and cl[0].length - cur_cluster[-1][0].length > max(MERGE_THR, cl[0].length * MERGE_RAT):
                            clusters.append(cur_cluster)
                            cur_cluster = [cl]
                        else:
                            cur_cluster.append(cl)
                    if cur_cluster:
                        clusters.append(cur_cluster)
                    for cl in clusters:
                        if len(cl) > 1:
                            new_len = int(np.mean([c[0].length for c in cl]))
                            by_genome_id = defaultdict(list)
                            bp_1 = cl[0][0].bp_1
                            bp_2 = cl[0][0].bp_2
                            if key == 'DEL':
                                bp_2.position = bp_1.position + new_len
                            for c in cl:
                                for db in c:
                                    by_genome_id[(db.genome_id, db.haplotype_1)].append(db)
                            for dbs in by_genome_id.values():
                                dbs[0].bp_1 = bp_1
                                dbs[0].bp_2 = bp_2
                                dbs[0].length = new_len
                                dbs[0].prec = 0
                                dbs[0].supp = sum([db.supp for db in dbs])
                                read_ids = []
                                for db in dbs:
                                    read_ids += db.supp_read_ids
                                dbs[0].supp_read_ids = list(set(read_ids))
                                for db in dbs[1:]:
                                    db.is_pass = 'FAIL_VNTR'
    
def get_phasingblocks(hb_vcf):
#    MIN_BLOCK_LEN = 10000
#    MIN_SNP = 10
    
    vcf = pysam.VariantFile(hb_vcf)
    haplotype_blocks = defaultdict(list)
    id_list = defaultdict(list)

    for var in vcf:
        if 'PS' in var.samples.items()[0][1].items()[-1] and var.samples.items()[0][1]['PS']:
            haplotype_blocks[(var.chrom, var.samples.items()[0][1]['PS'])].append(var.pos)

    phased_lengths = []
    for (chr_id, block_name), coords in haplotype_blocks.items():
        #if max(coords) - min(coords) > MIN_BLOCK_LEN and len(coords) >= MIN_SNP:
        phased_lengths.append(max(coords) - min(coords))
        id_list[chr_id].append(block_name)

    total_phased = sum(phased_lengths)
    _l50, n50 = _calc_nx(phased_lengths, total_phased, 0.50) 
    logger.info(f"\tTotal phased length: {total_phased}")
    logger.info(f"\tPhase blocks N50: {n50}")

    return id_list


def add_phaseset_id(double_breaks, id_list):
    for db in double_breaks:
        ind_1 = bisect.bisect_left(id_list[db.bp_1.ref_id], db.bp_1.position)
        ind_2 = bisect.bisect_left(id_list[db.bp_2.ref_id], db.bp_2.position)
        id1 = id_list[db.bp_1.ref_id][ind_1-1] if 0 < ind_1 < len(id_list[db.bp_1.ref_id]) else 0
        id2 = id_list[db.bp_2.ref_id][ind_2-1] if 0 < ind_2 < len(id_list[db.bp_2.ref_id]) else 0
        db.phaseset_id = (id1, id2)
    
def segment_coverage(histograms, genome_id, ref_id, ref_start, ref_end, haplotype):
    hist_start = ref_start // COV_WINDOW
    hist_end = ref_end // COV_WINDOW
    cov = [0,0,0]
    for i in [0,1,2]:
        cov_list = histograms[(genome_id, i, ref_id)][hist_start : hist_end + 1]
        if cov_list:
            cov[i] = int(np.median(cov_list))
    return (cov[haplotype],sum(cov))

def get_segments_coverage(db_segments, coverage_histograms, max_genomic_length):
    genomic_segments = defaultdict(list)
        
    for db, segments in db_segments.items():
        for (genome_id, seg_ref, seg_start, seg_end, seg_hp_ls, seg_hp) in segments:
            if seg_end - seg_start > max_genomic_length:
                if seg_start in [db.bp_1.position,  db.bp_2.position]:
                    seg_end = seg_start + max_genomic_length
                else: 
                    seg_start = seg_end - max_genomic_length
            
            (coverage, total_coverage) = segment_coverage(coverage_histograms, genome_id, seg_ref, seg_start, seg_end, seg_hp)
            
            if db.is_dup and seg_end - seg_start == 0:
                continue
                    
            gs = GenomicSegment(genome_id, seg_hp_ls, seg_ref, seg_start, seg_end,
                                coverage,total_coverage, seg_end - seg_start)
            genomic_segments[db].append(gs)
                
    return genomic_segments

def get_ref_adj(genomic_segments, ref_adj):
    adj_segments = []
    for key, pos_ls in ref_adj.items():
        if key[1] == 2:
            continue
        pos_ls2 = ref_adj[(key[0], 2)]
        for i, (pos1, db1) in enumerate(pos_ls):
            (pos2, db2) = pos_ls2[i]
            if db1 == db2:
                continue
            gs1 = genomic_segments[db1]
            gs2 = genomic_segments[db2]
            gs = [gs for gs in gs1 if gs.pos2 == pos1]
            g = [g for g in gs2 if g.pos1 == pos2]
            if gs and g:
                adj_segments.append((gs[0],g[0]))
    return adj_segments


def get_genomic_segments(double_breaks, coverage_histograms, hb_vcf, key_type, ref_lengths, min_ref_flank, max_genomic_length, min_sv_size):
    hb_points = []
    
    if hb_vcf:
        hb_points = get_phasingblocks(hb_vcf)
    
    if key_type == 'germline' and hb_points:
        add_phaseset_id(double_breaks, hb_points)
        cluster_inversions(double_breaks, coverage_histograms, min_sv_size)
        
    clusters = defaultdict(list)
    for br in double_breaks:
        clusters[br.to_string()].append(br)
        
    by_genome = defaultdict(list)
    for cl in clusters.values():
        gen_id = sorted(list(set([db.genome_id for db in cl])))
        gen_id = ','.join(gen_id)
        by_genome[gen_id] += cl
    db_segments = defaultdict(list)
    ref_adj = defaultdict(list)
    for d_breaks in by_genome.values():
        calc_gen_segments(d_breaks, coverage_histograms,ref_lengths, min_ref_flank, max_genomic_length, db_segments, ref_adj)
    genomic_segments = get_segments_coverage(db_segments, coverage_histograms, max_genomic_length)
    adj_segments = get_ref_adj(genomic_segments, ref_adj)
    return (genomic_segments, adj_segments)

def calc_gen_segments(double_breaks,coverage_histograms,ref_lengths, min_ref_flank, max_genomic_length, db_segments, ref_adj):

    bp1_list = defaultdict(list)
    bp2_list = defaultdict(list)
    single_bp = defaultdict(list)
    DEL_THR = 10000
    MAX_REF_DIST = 1000000
    BUFF = 50
    THR = 1
    for double_bp in double_breaks:
        if double_bp.bp_2.is_insertion or double_bp.sv_type == 'reciprocal_inv':
            continue
        if (double_bp.vcf_sv_type == 'DEL' or double_bp.vcf_sv_type == 'DUP') and double_bp.length < DEL_THR:
            continue
        db_cov = double_bp.supp + double_bp.bp_1.spanning_reads[double_bp.genome_id][double_bp.haplotype_1]
        bp1_list[(double_bp.genome_id, double_bp.bp_1.ref_id)].append(double_bp)
        single_bp[(double_bp.genome_id, double_bp.bp_1.ref_id)].append((double_bp.bp_1.position, double_bp.direction_1, '', (double_bp, double_bp.supp, double_bp.phaseset_id[0], double_bp.haplotype_1,db_cov)))
        if double_bp.is_dup:
            single_bp[(double_bp.genome_id, double_bp.bp_1.ref_id)].append((double_bp.bp_1.position, 1, 'dup', (double_bp, double_bp.supp, double_bp.phaseset_id[0], double_bp.haplotype_1,db_cov)))
        
        bp2_list[(double_bp.genome_id, double_bp.bp_2.ref_id)].append(double_bp)
        db_cov = double_bp.supp + double_bp.bp_2.spanning_reads[double_bp.genome_id][double_bp.haplotype_2]
        single_bp[(double_bp.genome_id, double_bp.bp_2.ref_id)].append((double_bp.bp_2.position, double_bp.direction_2, '', (double_bp, double_bp.supp, double_bp.phaseset_id[1], double_bp.haplotype_2,db_cov)))
        if double_bp.is_dup:
            single_bp[(double_bp.genome_id, double_bp.bp_2.ref_id)].append((double_bp.bp_2.position, -1, 'dup', (double_bp, double_bp.supp, double_bp.phaseset_id[1], double_bp.haplotype_2,db_cov)))
    for (genome_name, ref_name), s_bp in single_bp.items():
        for ps in s_bp:
            if ps[2]:
                ref_adj[ref_name,1].append((ps[0], ps[3][0]))
                ref_adj[ref_name,2].append((ps[0], ps[3][0]))
        s_bp = list(set(s_bp))
        s_bp.sort(key=lambda x: (x[0],x[1]))
        s_bp = [(min_ref_flank, 0, 0,(0,0,0,0,0))] + s_bp + [(ref_lengths[ref_name] - min_ref_flank, 0, 0,(0,0,0,0,0))]
        pos_ls = [s[0] for s in s_bp if not s[1] == -1]
        neg_ls = [s[0] for s in s_bp if not s[1] == 1]
        pos_bp = [s[3] for s in s_bp if not s[1] == -1]
        neg_bp = [s[3] for s in s_bp if not s[1] == 1]
        bp1 = bp1_list[(genome_name, ref_name)]
        
        for (seg1, seg2) in zip(s_bp[:-1], s_bp[1:]):
            if seg1[2] or seg2[2]:
                continue
            if seg1[1] == 1 and seg2[1] == -1 and seg2[0] - seg1[0] <= MAX_REF_DIST and (abs(seg1[3][1] - seg2[3][1]) <= min([seg1[3][1], seg2[3][1]])):
                if not (seg1[3][2] == seg2[3][2] and not seg1[3][3] == seg1[3][3] and 0 not in [seg1[3][3], seg1[3][3]]):
                    ref_adj[ref_name,1].append((seg1[0],seg1[3][0]))
                    ref_adj[ref_name,2].append((seg2[0],seg2[3][0]))
                
        for db in bp1:
            if db.is_single:
                continue
            bp1_len = np.median([c[0].segment_length for c in db.bp_1.connections if c in db.bp_2.connections])
            db_cov = db.supp + db.bp_1.spanning_reads[db.genome_id][db.haplotype_1]
            if db.bp_1.dir_1 == -1:
                ind = bisect.bisect_right(pos_ls,db.bp_1.position)
                if not db.is_dup and db.bp_1.position == pos_ls[ind]:
                    ind = ind +1
                ind_max = bisect.bisect_right(pos_ls, db.bp_1.position + max_genomic_length)
                ind_1 = [i for i in range(ind, min([ind_max + 2, len(pos_ls)])) if pos_bp[i][0] and set(db.supp_read_ids).intersection(set(pos_bp[i][0].supp_read_ids))]
                if not ind_1:
                    indb = bisect.bisect_right(pos_ls,db.bp_1.position + bp1_len - BUFF)
                    ind = min([max([ind, indb]), len(pos_ls)-1])
                    for i in range(ind, min([ind_max + 2, len(pos_ls)])):
                        if not pos_bp[i][0]:
                            ind_1.append(i)
                            continue
                        if not db.haplotype_1 == 0 and not pos_bp[i][3] == 0 and (db.phaseset_id[0] == pos_bp[i][2] and not db.haplotype_1 == pos_bp[i][3]):
                            continue
                        if not abs(db_cov - pos_bp[i][4]) < min(db_cov, pos_bp[i][4])*THR:
                            continue
                        ind_1.append(i)
                if ind_1:
                    pos2 = pos_ls[ind_1[0]]
                    hp2 = pos_bp[ind_1[0]][3]
                else:
                    pos2 = db.bp_1.position + max_genomic_length
                    hp2 = db.haplotype_1
                db_segments[db].append((genome_name, ref_name, db.bp_1.position, pos2, (db.haplotype_1, hp2),db.haplotype_1))
                
            elif db.bp_1.dir_1 == 1:
                ind = bisect.bisect_left(neg_ls,db.bp_1.position)
                if neg_ls[ind -1] == db.bp_1.position:
                    ind = ind -1
                ind_max = bisect.bisect_left(neg_ls , db.bp_1.position - max_genomic_length)
                ind_1 = [i for i in range(max(0, ind_max-1), ind) if neg_bp[i][0] and set(db.supp_read_ids).intersection(set(neg_bp[i][0].supp_read_ids))]
                if not ind_1:
                    indb = bisect.bisect_left(neg_ls,db.bp_1.position - bp1_len + BUFF)
                    ind = max([min([ind, indb]), 1 ])
                    for i in range(max(0, ind_max-1), ind):
                        if not neg_bp[i][0]:
                            ind_1.append(i)
                            continue
                        if not db.haplotype_1 == 0 and not neg_bp[i][3] == 0 and (db.phaseset_id[0] == neg_bp[i][2] and not db.haplotype_1 == neg_bp[i][3]):
                            continue
                        if not abs(db_cov - neg_bp[i][4]) < min(db_cov, neg_bp[i][4])*THR:
                            continue
                        ind_1.append(i)
                if ind_1:
                    pos2 = neg_ls[ind_1[-1]]
                    hp2 = neg_bp[ind_1[-1]][3]
                else:
                    pos2 = db.bp_1.position - max_genomic_length
                    hp2 = db.haplotype_1
                db_segments[db].append((genome_name, ref_name, pos2, db.bp_1.position, ( hp2,db.haplotype_1),db.haplotype_1))
                
        bp2 = bp2_list[(genome_name, ref_name)]
        for db in bp2:
            if db.is_single:
                continue
            bp2_len = np.median([c[1].segment_length for c in db.bp_1.connections if c in db.bp_2.connections])
            db_cov = db.supp + db.bp_2.spanning_reads[db.genome_id][db.haplotype_2]
            if db.bp_2.dir_1 == -1:
                ind = bisect.bisect_right(pos_ls,db.bp_2.position)
                if db.bp_2.position == pos_ls[ind]:
                    ind = ind +1
                ind_max = bisect.bisect_right(pos_ls, db.bp_2.position + max_genomic_length)
                ind_1 = [i for i in range(ind, min([ind_max + 2, len(pos_ls)])) if pos_bp[i][0] and set(db.supp_read_ids).intersection(set(pos_bp[i][0].supp_read_ids))]
                if not ind_1:
                    indb = bisect.bisect_right(pos_ls,db.bp_2.position + bp2_len)
                    ind = min([max([ind, indb]), len(pos_ls)-1])
                    for i in range(ind, min([ind_max + 2, len(pos_ls)])):
                        if not pos_bp[i][0]:
                            ind_1.append(i)
                            continue
                        if not db.haplotype_2 == 0 and not pos_bp[i][3] == 0 and (db.phaseset_id[1] == pos_bp[i][2] and not db.haplotype_2 == pos_bp[i][3]):
                            continue
                        if not abs(db_cov - pos_bp[i][4]) < min(db_cov, pos_bp[i][4])*THR:
                            continue
                        ind_1.append(i)
                if ind_1:
                    pos2 = pos_ls[ind_1[0]]
                    hp2 = pos_bp[ind_1[0]][3]
                else:
                    pos2 = db.bp_2.position + max_genomic_length
                    hp2 = db.haplotype_2
                db_segments[db].append((genome_name, ref_name, db.bp_2.position, pos2, (db.haplotype_2, hp2),db.haplotype_2))
                
            elif db.bp_2.dir_1 == 1:
                ind = bisect.bisect_left(neg_ls,db.bp_2.position)
                if not db.is_dup and neg_ls[ind -1] == db.bp_2.position:
                    ind = ind -1
                ind_max = bisect.bisect_left(neg_ls , db.bp_2.position - max_genomic_length)
                ind_1 = [i for i in range(max(0, ind_max-1), ind) if neg_bp[i][0] and set(db.supp_read_ids).intersection(set(neg_bp[i][0].supp_read_ids))]
                if not ind_1:
                    indb = bisect.bisect_left(neg_ls,db.bp_2.position - bp2_len)
                    ind = max([min([ind, indb]), 1 ])
                    for i in range(max(0, ind_max-1), ind):
                        if not neg_bp[i][0]:
                            ind_1.append(i)
                            continue
                        if not db.haplotype_2 == 0 and not neg_bp[i][3] == 0 and (db.phaseset_id[1] == neg_bp[i][2] and not db.haplotype_2 == neg_bp[i][3]):
                            continue
                        if not abs(db_cov - neg_bp[i][4]) < min(db_cov, neg_bp[i][4])*THR:
                            continue
                        ind_1.append(i)
                if ind_1:
                    pos2 = neg_ls[ind_1[-1]]
                    hp2 = neg_bp[ind_1[-1]][3]
                else:
                    pos2 = db.bp_2.position - max_genomic_length
                    hp2 = db.haplotype_2
                db_segments[db].append((genome_name, ref_name, pos2, db.bp_2.position, (hp2, db.haplotype_2), db.haplotype_2))        

    
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

def match_haplotypes(double_breaks):
    PHASE_THR = 0.66
    
    clusters = defaultdict(list)
    for br in double_breaks:
        if br.is_pass == 'PASS':
            clusters[br.to_string()].append(br)
            
    for cl in clusters.values():
        by_genome_id = defaultdict(list)
        for db in cl:
            by_genome_id[db.genome_id].append(db)
            
        for genome_id, db_ls in by_genome_id.items():
            if len(db_ls) == 1:
                continue
            
            by_haplotype1 = defaultdict(list)
            haplotype1_supp = defaultdict(int)
            for db in db_ls:
                by_haplotype1[db.haplotype_1].append(db)
                haplotype1_supp[db.haplotype_1] += db.supp
            hp_list1 = sorted(haplotype1_supp, key=lambda k:haplotype1_supp[k], reverse=True)
            hp_list1_phased = [x for x in hp_list1 if not x == 0]
            
            if 0 in hp_list1 and hp_list1_phased:
                by_haplotype1[hp_list1_phased[0]] += by_haplotype1[0]
                by_haplotype1[0] = []
                
            if len(hp_list1_phased) > 1:
                if haplotype1_supp[hp_list1_phased[0]] * PHASE_THR > haplotype1_supp[hp_list1_phased[1]]:
                    by_haplotype1[hp_list1_phased[0]] += by_haplotype1[hp_list1_phased[1]]
                    by_haplotype1[hp_list1_phased[1]] = []
                    
            for hp1, db_list2 in by_haplotype1.items():
                if not db_list2:
                    continue
                
                by_haplotype2 = defaultdict(list)
                haplotype2_supp = defaultdict(int)
                for db in db_list2:
                    by_haplotype2[db.haplotype_2].append(db)
                    haplotype2_supp[db.haplotype_2] += db.supp
                    
                hp_list2 = sorted(haplotype2_supp, key=lambda k:haplotype2_supp[k], reverse=True)
                hp_list2_phased = [x for x in hp_list2 if not x == 0]
                
                if 0 in hp_list2 and hp_list2_phased:
                    by_haplotype2[hp_list2_phased[0]] += by_haplotype2[0]
                    by_haplotype2[0] = []
                    
                if len(hp_list2_phased) > 1:
                    if haplotype2_supp[hp_list2_phased[0]] * PHASE_THR > haplotype2_supp[hp_list2_phased[1]]:
                        by_haplotype2[hp_list2_phased[0]] += by_haplotype2[hp_list2_phased[1]]
                        by_haplotype2[hp_list2_phased[1]] = []
                        
                for hp2, dbs in by_haplotype2.items():
                    if not dbs:
                        continue
                    sum_supp = []
                    for db in dbs:
                        sum_supp += db.supp_read_ids
                    hp1_list = [db.haplotype_1 for db in dbs]
                    hp2_list = [db.haplotype_2 for db in dbs]
                    db = dbs[0]
                    db.supp_read_ids = list(set(sum_supp))
                    db.supp = len(db.supp_read_ids)
                    db.haplotype_1 = hp1
                    db.haplotype_2 = hp2
                    db.haplotypes = [hp1_list, hp2_list]
                    for db in dbs[1:]:
                        db.is_pass = 'FAIL_MERGED_HP'
                        
def cluster_db(db_list, coverage_histograms, min_sv_size):
    clusters = defaultdict(list)
    for br in db_list:
        br.subgraph_id = []
        if br.bp_1.is_insertion or not br.is_pass == 'PASS':
            continue
        clusters[br.to_string()].append(br)
    
    by_genome = defaultdict(list)
    for cl in clusters.values():
        gen_id = sorted(list(set([db.genome_id for db in cl])))
        gen_id = ','.join(gen_id)
        by_genome[gen_id] += cl
    
    for ind_id, db_list in enumerate(by_genome.values()):
        clusters = defaultdict(list) 
        for br in db_list:
            clusters[br.to_string()].append(br)
        conn_duplications(clusters, coverage_histograms)
        conn_inter(clusters, ind_id)
         
def cluster_inversions(double_breaks, coverage_histograms, min_sv_size):
    clusters = defaultdict(list) 
    for br in double_breaks:
        if br.bp_1.dir_1 == br.bp_2.dir_1 and not br.is_single:
            clusters[br.to_string()].append(br)
        
    by_genome = defaultdict(list)
    for cl in clusters.values():
        gen_id = sorted(list(set([db.genome_id for db in cl])))
        gen_id = ','.join(gen_id)
        by_genome[gen_id] += cl
    
    for ind_id, db_list in enumerate(by_genome.values()):
        complex_inv(db_list, coverage_histograms, min_sv_size, ind_id)

def reciprocal_inv(clusters):
    inv_list = defaultdict(list)
    MINSIZE = 200
    for cl in clusters.values():
        db = cl[0]
        supp = sum([db.supp for db in cl])
        pos_list = (cl[0].bp_1.position, cl[0].bp_1.dir_1, cl[0].haplotype_1, cl[0].phaseset_id[0], cl, supp)
        inv_list[db.bp_1.ref_id].append(pos_list)
        if not db.bp_1.ref_id == db.bp_2.ref_id:
            pos_list = (cl[0].bp_2.position, cl[0].bp_2.dir_1, cl[0].haplotype_2, cl[0].phaseset_id[1], cl, supp)
            inv_list[db.bp_2.ref_id].append(pos_list)
            
    for seq, pos_list in inv_list.items():
        pos_list.sort(key=lambda s:(s[0], -s[1]))
        neg_ls = [s for s in pos_list if s[1] == -1 and s[4][0].bp_1.ref_id == s[4][0].bp_2.ref_id]
        neg_pos = [s[0] for s in neg_ls]
        for i, pos in enumerate(pos_list[:-1]):
            if pos[1] == -1 or not pos[4][0].bp_1.ref_id == pos[4][0].bp_2.ref_id:
                continue
            strt = bisect.bisect_left(neg_pos, pos[0] - MINSIZE)
            end = min(bisect.bisect_left(neg_pos, pos[0] + MINSIZE), len(neg_pos)-2)
            black_list = []
            if len(pos[4]) == 1 and not pos_list[i][2] == 0:
                black_list = [ind for ind, (_,dir1, hp, phase_id, _, _) in enumerate(neg_ls[strt:end+1]) if phase_id == pos[3] and not (hp == 0 or hp == pos[2])]
            for j in range(strt,end+1):
                if not j in black_list and neg_ls[j][4][0].bp_2.position + MINSIZE >= pos[4][0].bp_2.position:
                    if (pos[4][0].bp_2.position - neg_ls[j][4][0].bp_1.position) / (neg_ls[j][4][0].bp_2.position - pos[4][0].bp_1.position) > 0.90:
                        for db in pos[4]+neg_ls[j][4]:
                            db.sv_type = 'reciprocal_inv'
    
    
def complex_inv(double_breaks, coverage_histograms, min_sv_size, ind_id):
    THR_MIN = 4
    THR_MAX = 10
    inv_list = defaultdict(list)
    MAX_DIST = 2000000
    
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
    
    reciprocal_inv(clusters)
    foldback_inv(clusters, coverage_histograms, ind_id)
    
    for cl in clusters.values():
        db = cl[0]
        if db.sv_type:
            continue
        
        supp = sum([db.supp for db in cl])
        pos_list = (cl[0].bp_1.position, cl[0].bp_1.dir_1, cl[0].haplotype_1, cl[0].phaseset_id[0], cl, supp)
        inv_list[db.bp_1.ref_id].append(pos_list)
        if not db.bp_1.ref_id == db.bp_2.ref_id:
            pos_list = (cl[0].bp_2.position, cl[0].bp_2.dir_1, cl[0].haplotype_2, cl[0].phaseset_id[1], cl, supp)
            inv_list[db.bp_2.ref_id].append(pos_list)
    
    for seq, pos_list in inv_list.items():

        pos_list.sort(key=lambda s:(s[0], -s[1]))
        pairs = []
        pairs_pos = []
        for i, pos in enumerate(pos_list[:-1]):
            supp = pos[5]
            supp_thr = min(max(supp * 0.5, THR_MIN) ,THR_MAX)
            black_list = [i]
            j_lis = []
            max_dist = pos[0] + MAX_DIST
            if len(pos[4]) == 1 and not pos_list[i][2] == 0:
                black_list += [ind for ind, (_,dir1, hp, phase_id, _, _) in enumerate(pos_list) if phase_id == pos[3] and not (hp == 0 or hp == pos[2])]
            for j in range (i+1, len(pos_list)):
                if not j in black_list and supp - supp_thr <= pos_list[j][5] <= supp + supp_thr and pos_list[j][1] + pos[1] == 0 and pos_list[j][0] <= max_dist:
                    j_lis.append(j)
            if not j_lis:
                continue
            if len(j_lis) == 1:
                jk = j_lis[0]
            elif len(j_lis) > 1:
                db1 = pos[4][0]
                pos1 = (db1.bp_1.position, db1.bp_2.position) if not db1.bp_1.ref_id == db1.bp_2.ref_id else (db1.bp_1.position, db1.bp_1.position)
                dbs = [pos_list[j][4][0] for j in j_lis]
                dist = []
                for db in dbs:
                    pos2 = (db.bp_1.position, db.bp_2.position) if not db.bp_1.ref_id == db.bp_2.ref_id else (db.bp_1.position, db.bp_1.position)
                    dist.append(min([abs(pos1[0] - pos2[0]), abs(pos1[1] - pos2[1])]))
                jk = j_lis[np.argmin(dist)]
            pairs.append((pos[4],pos_list[jk][4]))
            pairs_pos.append((pos[0],pos_list[jk][0]))
        
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
                        sv_type = 'Templated_ins_inv'
            
            for cl in cl_list:
                for db in cl:
                    #db.gr_id = 'INV' + str(ind_id)+ '_' + seq + '_' + str(i)
                    db.sv_type = sv_type
            
                
def foldback_inv(clusters, coverage_histograms, ind_id):

    inv_list_pos = defaultdict(list)
    inv_list_neg = defaultdict(list)
    FOLD_BACK_THR = 0.5
    FOLD_BACK_DIST_THR = 50000
    MINSIZE = 100000
    t = 0
    HP = [0,1,2]
    for cl in clusters.values():
        db = cl[0]
        if db.sv_type:
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
            pos1 = db.bp_1.position // COV_WINDOW
            cov1_0 = sum([coverage_histograms[(db.genome_id, hp, db.bp_1.ref_id)][pos1-1] for hp in HP])
            cov1_2 = sum(db.bp_1.spanning_reads[db.genome_id])
            cov2 = sum(db.bp_2.spanning_reads[db.genome_id])
            thr = int(db.supp * FOLD_BACK_THR)
            if cov1_0 > cov1_2 + thr and  cov1_2 > cov2 + thr:
                foldback_pairs['pos'].append(ind)
                for db in pos_lis[ind]:
                    db.sv_type = 'foldback'
        for ind, cl in enumerate(neg_lis):
            db = cl[0]
            if db.bp_2.position - db.bp_1.position > FOLD_BACK_DIST_THR:
                continue
            pos1 = db.bp_1.position // COV_WINDOW
            cov1_0 = sum([coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos1-1] for hp in HP])
            cov1_2 = sum(db.bp_1.spanning_reads[db.genome_id])
            cov2 = sum(db.bp_2.spanning_reads[db.genome_id])
            thr = int(db.supp * FOLD_BACK_THR)
            if cov2 > cov1_2 + thr and  cov1_2 > cov1_0 + thr:
                foldback_pairs['neg'].append(ind)
                for db in neg_lis[ind]:
                    db.sv_type = 'foldback'
        if len(foldback_pairs) == 2:
            pos_ls = []
            neg_ls = []
            all_ls = [db for ind in foldback_pairs['pos'] for db in pos_lis[ind]] + [db for ind in foldback_pairs['neg'] for db in neg_lis[ind]]
            all_ls.sort(key=lambda b:b.bp_1.position)
            neg = True
            for i, db in enumerate(all_ls):
                if neg:
                    if db.bp_1.dir_1 == -1:
                        neg_ls.append(db)
                    elif i == len(all_ls) -1 or all_ls[i+1].bp_1.dir_1 == 1:
                        neg = False
                        pos_ls.append(db)
                else:
                    if db.bp_1.dir_1 == 1:
                        pos_ls.append(db)
            if pos_ls and neg_ls and pos_ls[-1].bp_2.position - neg_ls[0].bp_1.position > MINSIZE:
                supp_pos = [db.supp for db in pos_ls]
                supp_neg = [db.supp for db in neg_ls]
                if abs(sum(supp_pos)- sum(supp_neg)) <= min(supp_pos + supp_neg):
                    pos1 = neg_ls[0].bp_1.position// COV_WINDOW
                    pos2 = pos_ls[-1].bp_2.position//COV_WINDOW
                    seg_cov = int(np.median([sum([coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos] for hp in HP]) for pos in range(pos1,pos2)]))
                    seg2 = int(np.median([sum([coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos] for hp in HP]) for pos in range(pos2,pos2+5)]))
                    seg1 = int(np.median([sum([coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos] for hp in HP]) for pos in range(pos1-5,pos1)]))
                    if seg_cov > max(seg1,seg2) * 1.5 and abs(seg1 - seg2) >= min(seg1,seg2)*0.75:
                        t+=1
                        for db in all_ls:
                            db.sv_type = 'BFB_foldback'
                            db.gr_id = 'BFB' +  str(ind_id)+ '_' + str(t)                        

def cluster_indels(double_breaks):
    by_genome = defaultdict(list)
    DEL_THR = 10000
    components_list = []
    for db in double_breaks:
        if db.vcf_sv_type == 'DEL' or db.vcf_sv_type =='DUP' and db.length < DEL_THR:
            by_genome[(db.genome_id, db.mut_type, db.bp_1.ref_id)].append(db)
    for db_list in by_genome.values():
        components_list += indel_clust(db_list)
        
    return components_list
    
  
def indel_clust(db_list):
    DIFF_THR = 20000
    SIZE_THR = 5
    LEN_THR = 2
    cur_cluster = []
    clusters = []
    components_list = []
    junction_type = defaultdict(int)
    for db in db_list:
        if cur_cluster and db.bp_1.position - cur_cluster[-1].bp_1.position > DIFF_THR:
            if len(cur_cluster) > SIZE_THR:
                clusters.append(cur_cluster)
            cur_cluster = [db]
        else:
            cur_cluster.append(db)
    if cur_cluster:
        clusters.append(cur_cluster)
    for cl in clusters:
        t = 0
        sum_len = 0
        for db1, db2 in zip(cl[:-1], cl[1:]):
            if db2.bp_1.position - db1.bp_1.position > db1.length:
                t+=1
                sum_len +=db1.length
        if t > SIZE_THR and sum_len * LEN_THR > (cl[-1].bp_1.position - cl[0].bp_1.position):
            sv_type = defaultdict(int)
            for db in cl:
                sv_type[db.vcf_sv_type]+=1
                jtype = 'TH' if str(db.direction_1) + str(db.direction_2) == '-11' else 'HT'
                junction_type[jtype] +=1
            components_list.append((sv_type, junction_type,'indel', 1, cl, cl[0].genome_id))
    return components_list
            
def conn_duplications(clusters, coverage_histograms):
    
    DUP_COV_THR = 0.5
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
            cov3 = int(np.median(coverage_histograms[(db.genome_id, db.haplotype_1, db.bp_1.ref_id)][pos1:min(pos2+1, max_pos)]))
            if cov3 > cov1 + db.supp * DUP_COV_THR:
                is_dup = True
            svtype = 'tandem_duplication' if is_dup else 'Templated_ins'
            for db in cl:
                db.is_dup = is_dup
                db.sv_type = svtype
            
def conn_inter(clusters, ind_id):
    intra_chr_cl = defaultdict(list)
    DEL_THR = 1000000
    t = 0
    THR_MIN = 4
    THR_MAX = 10
    BUFF = 100
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
            supp = sum([db.supp for db in cl])
            supp_thr = min(max(THR_MIN, supp * 0.5), THR_MAX)
            dir1, dir2 = -1 * db.direction_1, -1 * db.direction_2
            for ind2, cl2 in enumerate(intra_chr[ind+1:]):
                db2 = cl2[0]
                supp2 = sum([db2.supp for db2 in cl2])
                sv_type = ''
                if (db2.direction_1, db2.direction_2) == (dir1, dir2) and abs(supp - supp2) <= supp_thr:
                    dist = max([abs(db.bp_1.position - db2.bp_1.position), abs(db.bp_2.position - db2.bp_2.position)]) 
                    if dist > DEL_THR:
                        continue
                    
                    bp1_conn = [c for c in db.supp_read_ids if c in db2.supp_read_ids]
                    if bp1_conn:
                        sv_type = 'Templated_ins'
                    else:
                        chr_ls = defaultdict(list)
                        chr_ls[db.bp_1.ref_id].append(db.bp_1.dir_1 * db.bp_1.position)
                        chr_ls[db.bp_2.ref_id].append(db.bp_2.dir_1 * db.bp_2.position)
                        chr_ls[db2.bp_1.ref_id].append(db2.bp_1.dir_1 * db2.bp_1.position)
                        chr_ls[db2.bp_2.ref_id].append(db2.bp_2.dir_1 * db2.bp_2.position)
                        
                        seglen_ls = defaultdict(list)
                        seglen_ls[db.bp_1.ref_id].append(max([c[0].segment_length for c in db.bp_1.connections if c in db.bp_2.connections]))
                        seglen_ls[db.bp_2.ref_id].append(max([c[1].segment_length for c in db.bp_1.connections if c in db.bp_2.connections]))
                        seglen_ls[db2.bp_1.ref_id].append(max([c[0].segment_length for c in db2.bp_1.connections if c in db2.bp_2.connections]))
                        seglen_ls[db2.bp_2.ref_id].append(max([c[1].segment_length for c in db2.bp_1.connections if c in db2.bp_2.connections]))
                    
                        rec = 0
                        for key, seglen in chr_ls.items():
                            if sum(seglen) <= 0:
                                rec += 1
                                continue
                            if max(seglen_ls[key]) > sum(seglen) + BUFF:
                                rec +=1
                        if rec == 2:
                            sv_type = 'Reciprocal_tra'
                        else:
                            sv_type = 'Templated_ins'
                
                        ind_list.append(ind2)
                        t += 1
                        for db in cl:
                            db.sv_type = sv_type
                            db.gr_id = 'TIC' + str(ind_id) + '_' + str(t)
                            
                        for db2 in cl2:
                            db2.sv_type = sv_type
                            db2.gr_id = 'TIC' + str(ind_id) + '_' + str(t)
                if sv_type:
                    continue

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
        
    logger.info('Extracting split alignments')
    split_reads = get_splitreads(segments_by_read)
    ins_list_all = get_insertionreads(segments_by_read)
    cont_id  = list(control_id)[0] if control_id else '' 
    
    logger.info('Extracting clipped reads')
    clipped_clusters = []
    extract_clipped_end(segments_by_read)
    clipped_reads = get_clipped_reads(segments_by_read)
    clipped_clusters = cluster_clipped_ends(clipped_reads, args.bp_cluster_size,args.min_ref_flank, ref_lengths)
    
    logger.info('Starting breakpoint detection')
    double_breaks, single_bps = get_breakpoints(split_reads, ref_lengths, args)
    
    logger.info('Clustering unmapped insertions')
    ins_clusters = extract_insertions(ins_list_all, clipped_clusters, ref_lengths, args)

    match_long_ins(ins_clusters, double_breaks, args.min_sv_size, args.tra_to_ins)
    
    if args.single_bp:
        logger.info('Starting single breakpoint detection')
        single_bps = get_single_bp(single_bps, clipped_clusters, double_breaks+ins_clusters, args.bp_min_support, cont_id,args.min_ref_flank, ref_lengths)
    else:
        single_bps = []
    
    logger.info('Starting compute_bp_coverage')
    if args.vntr_file:
        add_vntr_annot(double_breaks + ins_clusters, args)
    get_coverage_parallel(bam_files, genome_ids, thread_pool, args.min_mapping_quality, double_breaks + ins_clusters + single_bps)

        
    logger.info('Filtering breakpoints')
    double_breaks = double_breaks_filter(double_breaks, single_bps, args.bp_min_support, cont_id, args.resolve_overlaps, args.sv_size)
    double_breaks.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position, b.direction_1))
    
    if args.single_bp and single_bps:
        single_bps = filter_single_bp(single_bps, cont_id, args.control_vaf, args.vaf_thr)
    insertion_filter(ins_clusters, args.bp_min_support, cont_id)
    ins_clusters.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position))
   
    double_breaks +=  ins_clusters
    annotate_mut_type(double_breaks, cont_id, args.control_vaf, args.vaf_thr)
        
    logger.info('Writing breakpoints')
    output_breaks(double_breaks, genome_ids, args.phase_vcf, open(os.path.join(args.out_dir,"breakpoints_double.csv"), "w"))
    
    double_breaks = filter_fail_double_db(double_breaks, single_bps, coverage_histograms, segments_by_read, bam_files, thread_pool, args)
    
    
    if args.output_read_ids:
            output_readids(double_breaks['germline'], genome_ids, open(os.path.join(args.out_dir,"read_ids.csv"), "w"))
    
    return double_breaks
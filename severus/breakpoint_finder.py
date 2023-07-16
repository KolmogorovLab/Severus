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

from severus.bam_processing import _calc_nx, extract_clipped_end


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
                 "genome_id", 'bp_list', 'is_pass1', 'is_pass2', 'mapq_1', 'mapq_2', 'is_dup', 'has_ins')
    def __init__(self, ref_id_1, pos_1, pos_1_ori , sign_1, ref_id_2, pos_2, pos_2_ori, sign_2, haplotype_1, 
                 haplotype_2, read_id, genome_id, is_pass1, is_pass2, mapq_1, mapq_2, is_dup, has_ins):
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
                 'DR', 'DV', 'hvaf', 'vaf', 'prec')
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
    def to_string(self):
        strand_1 = "+" if self.direction_1 > 0 else "-"
        strand_2 = "+" if self.direction_2 > 0 else "-"
        label_1 = "{0}{1}:{2}".format(strand_1, self.bp_1.ref_id, self.bp_1.position)
        if self.bp_2.is_insertion:
            label_2 = "{0}:{1}".format('INS', self.length)
        elif self.bp_2.contig_id or self.bp_2.loose_end_id:
            label_2 = "{0}:{1}".format(self.bp_2.ref_id, self.bp_2.position)
        else:
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

    def left_coord_str(self):
        return f"{self.ref_id}:{self.pos1}"

    def right_coord_str(self):
        return f"{self.ref_id}:{self.pos2}"

    def left_coord_tuple(self):
        return (self.ref_id, self.pos1, self.dir1)

    def right_coord_tuple(self):
        return (self.ref_id, self.pos2, self.dir2)


def get_breakpoints(split_reads, ref_lengths, args):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """
    clust_len = args.bp_cluster_size
    min_reads = args.bp_min_support
    min_ref_flank = args.min_ref_flank 
    sv_size = args.min_sv_size
    MAX_SEGMENT_DIST= 500
    
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
        dist = abs(s2.read_start - s1.read_end)
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
                                s1.is_pass, s2.mapq, s1.mapq, is_dup, has_ins)
            seq_breakpoints_r[s2.ref_id].append(rc)
            seq_breakpoints_l[s1.ref_id].append(rc)
        else:
            if sign_1 == -1 and sign_2 == 1 and s2.ref_start <= s1.ref_start <= s2.ref_end <= s1.ref_end:
                is_dup = True
            rc = ReadConnection(s1.ref_id, ref_bp_1, ref_bp_1_ori, sign_1, s2.ref_id, ref_bp_2, ref_bp_2_ori, sign_2,
                                s1.haplotype, s2.haplotype, s1.read_id, s1.genome_id, s1.is_pass, 
                                s2.is_pass, s1.mapq, s1.mapq, is_dup, has_ins)
            seq_breakpoints_r[s1.ref_id].append(rc)
            seq_breakpoints_l[s2.ref_id].append(rc)
        
    for read_segments in split_reads:
        read_segments.sort(key=lambda x:x.read_start)
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            if abs(s2.read_start - s1.read_end) < MAX_SEGMENT_DIST:
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
            
        if max(by_genome_id.values()) >= min_reads:
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
    
    for x in cl:
        unique_reads[(x.genome_id,x.haplotype_1,x.haplotype_2)].add(x.read_id)
        if x.is_pass1 == 'PASS' and x.is_pass2 == 'PASS':
            unique_reads_pass[(x.genome_id,x.haplotype_1,x.haplotype_2)].add(x.read_id)
            
    by_genome_id_pass = defaultdict(int)
    happ_support_1 = defaultdict(list)
    happ_support_2 = defaultdict(list)
    for key, values in unique_reads.items():
        if unique_reads_pass[key]:
            by_genome_id_pass[key[0]] += len(unique_reads_pass[key])
            happ_support_1[key[0]].append(key[1])
            happ_support_2[key[0]].append(key[2])
            
    if by_genome_id_pass.values() and max(by_genome_id_pass.values()) >= min_reads:
        for keys in unique_reads.keys():
            genome_id = keys[0]
            haplotype_1 = keys[1]
            haplotype_2 = keys[2]
            
            if sum(happ_support_1[genome_id]) == 3 or sum(happ_support_2[genome_id]) == 3:
                genotype = 'hom'
            else:
                genotype = 'het'
                
            supp = len(unique_reads[keys])
            support_reads = unique_reads[keys]
            
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
            
    return db_list

#TODO BAM/HAPLOTYPE SPECIFIC FILTER
def double_breaks_filter(double_breaks, min_reads, genome_ids):

    PASS_2_FAIL_RAT = 0.5
    CONN_2_PASS = 0.5
    CHR_CONN = 2
    COV_THR = 3
    NUM_HAPLOTYPES = 3
    
    for db in double_breaks:
        conn_1 = [cn for cn in db.bp_1.connections if cn.genome_id == db.genome_id and cn.haplotype_1 == db.haplotype_1]
        conn_2 = [cn for cn in db.bp_2.connections if cn.genome_id == db.genome_id and cn.haplotype_2 == db.haplotype_2]#
        
        conn_pass_1 =[cn for cn in conn_1 if cn.is_pass1 == 'PASS']
        conn_pass_2 =[cn for cn in conn_2 if cn.is_pass2 == 'PASS']#
        
        conn_count_1 = Counter([cn.is_pass1 for cn in conn_1])
        conn_count_2 = Counter([cn.is_pass2 for cn in conn_2])#
        
        is_dup = [cn.is_dup for cn in conn_pass_1 if cn in conn_pass_2]
        
        if any(is_dup):
            db.is_dup = 'DUP'
            
        has_ins = [cn.has_ins for cn in conn_pass_1 if cn in conn_pass_2]
        if any(has_ins):
            db.has_ins = np.median(has_ins)
        
        if not conn_count_1['PASS'] or not conn_count_2['PASS']:
            db.is_pass = 'FAIL_LOWSUPP'
            continue#
            
        if conn_count_1['PASS'] < len(conn_1) * PASS_2_FAIL_RAT or conn_count_2['PASS'] < len(conn_2) * PASS_2_FAIL_RAT:
            db.is_pass = 'FAIL_MAP_CONS'
            continue#
            
        conn_valid_1 = Counter([len(cn.bp_list) for cn in conn_pass_1])
        conn_valid_2 = Counter([len(cn.bp_list) for cn in conn_pass_2])
        if conn_valid_1[2] < len(conn_pass_1) * CONN_2_PASS and conn_valid_2[2] < len(conn_pass_2) * CONN_2_PASS:
            db.is_pass = 'FAIL_CONN_CONS'
            continue#
            
        conn_ref_1 = Counter([cn.ref_id_1 for cn in conn_pass_1])
        conn_ref_2 = Counter([cn.ref_id_2 for cn in conn_pass_2])
        if len(conn_ref_1) > CHR_CONN or len(conn_ref_2) > CHR_CONN :
            db.is_pass = 'FAIL_CONN_CONS'
            continue#
            
        if db.supp < min_reads:
            db.is_pass = 'FAIL_LOWSUPP'
            continue#
            
    db_list = []
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
        
    for cl in clusters.values():
        count_pass = Counter([db1.is_pass for db1 in cl])
        if not count_pass['PASS']:
            continue
        
        for db in cl:
            db.is_pass = 'PASS'
        
        gen_ids = list(set(genome_ids) - set([db1.genome_id for db1 in cl]))
        if gen_ids:
            span_bp1 = defaultdict(int)
            span_bp2 = defaultdict(int)
            for gen_id in gen_ids:
                for i in range(NUM_HAPLOTYPES):
                    span_bp1[gen_id] += cl[0].bp_1.spanning_reads[(gen_id, i)]
                    span_bp2[gen_id] += cl[0].bp_2.spanning_reads[(gen_id, i)]
                
            for (genome_id, haplotype), count in cl[0].bp_1.spanning_reads.items():
                count = count if not haplotype == 0 else span_bp1[genome_id] 
                if genome_id in gen_ids and count < COV_THR:
                    for db in cl:
                        if db.haplotype_1 == haplotype:
                            db.is_pass = 'FAIL_LOWCOV_OTHER'
                            
            for (genome_id, haplotype), count in cl[0].bp_2.spanning_reads.items():
                count = count if not haplotype == 0 else span_bp2[genome_id] 
                if genome_id in gen_ids and count < COV_THR:
                    for db in cl:
                        if db.haplotype_2 == haplotype:
                            db.is_pass = 'FAIL_LOWCOV_OTHER'
                        
        db_list += cl
        
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
                    if cl_ins and abs(cl1.segment_length - cl_ins[-1].segment_length) > sv_len_diff:
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
            prec = 1
            for x in cl:
                unique_reads[(x.genome_id, x.haplotype)].add(x)
                if x.is_pass == 'PASS':
                    unique_reads_pass[(x.genome_id, x.haplotype)].add(x)#
                    
            by_genome_id_pass = defaultdict(int)
            happ_support_1 = defaultdict(list)#
            
            for key, values in unique_reads.items():
                if unique_reads_pass[key]:
                    by_genome_id_pass[key[0]] = len(set([red.read_id for red in unique_reads_pass[key]]))
                    happ_support_1[key[0]].append(key[1])#
                    
            if not by_genome_id_pass.values() or max(by_genome_id_pass.values()) < MIN_FULL_READ_SUPP:
                continue
            
            pos_list = [x.ref_end_ori for x in cl if x.is_pass == 'PASS']
            s_len = [x.segment_length for x in cl if x.is_pass == 'PASS']
            if np.median(s_len) * (len(s_len) + 1) < pos_list[-1] - pos_list[0]:
                continue
            
            if max(pos_list) - min(pos_list) > 2 * sv_size:
                prec = 0
            elif np.std(s_len) / np.mean(s_len) > CV_THR:
                prec = 0
                
            position = int(np.median(pos_list))
            mapq = int(np.median([x.mapq for x in cl if x.is_pass == 'PASS']))
            ins_length = int(np.median(s_len))
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
                    unique_reads = unique_reads_pass[key]
                    bp_1 = Breakpoint(seq, position, -1, mapq)
                    bp_1.read_ids = [x.read_id for x in unique_reads]
                    bp_1.connections = unique_reads
                    bp_3 = Breakpoint(seq, position, 1, mapq)
                    bp_3.is_insertion = True
                    bp_3.insertion_size = ins_length
                    supp = len(unique_reads)
                    supp_reads = unique_reads
                    genome_id = key[0]
                    if sum(happ_support_1[genome_id]) == NUM_HAPLOTYPES:
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


def insertion_filter(ins_clusters, min_reads, genome_ids):
    PASS_2_FAIL_RAT = 0.9
    COV_THR = 3
    NUM_HAPLOTYPES = 3
    
    for ins in ins_clusters:
        conn_1 = [cn for cn in ins.bp_1.connections if cn.genome_id == ins.genome_id and cn.haplotype == ins.haplotype_1]
        conn_count_1 = Counter([cn.is_pass for cn in conn_1])
        
        if not conn_count_1['PASS']:
            ins.is_pass = 'FAIL_MAP_CONS'
            continue
        
        if conn_count_1['PASS'] < len(conn_1) * PASS_2_FAIL_RAT:
            ins.is_pass = 'FAIL_MAP_CONS'
            continue
        
        if ins.supp < min_reads:
            ins.is_pass = 'FAIL_LOWSUPP'
            continue
        
    cur_cluster = []
    clusters = []
    ins_list = []
    for ins in ins_clusters:
        if cur_cluster and ins.bp_1.position == cur_cluster[-1].bp_1.position and ins.length == cur_cluster[-1].length:
            cur_cluster.append(ins)
        else:
            clusters.append(cur_cluster)
            cur_cluster = [ins]
    if cur_cluster:
        clusters.append(cur_cluster)#
    clusters = clusters[1:] 
        
    for cl in clusters:
        count_pass = Counter([db1.is_pass for db1 in cl])
        if not count_pass['PASS']:
            continue
        gen_ids = list(set(genome_ids) - set([db1.genome_id for db1 in cl]))
        if gen_ids:
            span_bp1 = defaultdict(int)
            for gen_id in gen_ids:
                for i in range(NUM_HAPLOTYPES):
                    span_bp1[gen_id] += cl[0].bp_1.spanning_reads[(gen_id, i)]
            for (genome_id, haplotype), count in cl[0].bp_1.spanning_reads.items():
                count = count if not haplotype == 0 else span_bp1[genome_id]
                if genome_id in gen_ids and count < COV_THR:
                    for ins in cl:
                        if ins.haplotype_1 == haplotype:
                            ins.is_pass = 'FAIL_LOWCOV_OTHER'
                        
        for ins in cl:
            ins_list.append(ins)
            
    return ins_list

def get_clipped_reads(segments_by_read):
    clipped_reads = defaultdict(list)
    for read in segments_by_read.values():
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

def tra_to_ins(ins_list_pos, ins_list, bp, dir_bp, dbs, ins_clusters, double_breaks, min_sv_size):
    
    INS_WIN = 2000 
    NUM_HAPLOTYPE = 3
    
    ins_1 = ins_list_pos[bp.ref_id]
    strt = bisect.bisect_left(ins_1, bp.position - INS_WIN)
    end = bisect.bisect_left(ins_1, bp.position + INS_WIN)
    ins_to_remove = []
    
    if strt == end:
        return []
    
    ins_db = ins_list[bp.ref_id]
    clusters = defaultdict(list)
    for ins in ins_db[strt:end]:
        clusters[ins.to_string()].append(ins)
    
    for ins_cl in clusters.values():
        gen_id_1 = defaultdict(list)
        hp_list = defaultdict(list)
        ins = ins_cl[0]
        if (ins.length + min_sv_size) < abs(ins.bp_1.position - bp.position):
            continue
        
        if dbs[0].bp_1.ref_id == dbs[0].bp_2.ref_id:
            svtype = 'Intra_chr_ins'
        else:
            svtype = 'Inter_chr_ins'
            
        for ins in ins_cl:
            hp_list[ins.genome_id].append(ins.haplotype_1)
        
        for db in dbs:
            gen_id_1[(db.genome_id, db.haplotype_1)].append(db)
            hp_list[db.genome_id].append(db.haplotype_1)
            db.sv_type = svtype
            
        for ins in ins_cl:
            ins_to_remove.append(ins)
            genotype = 'hom' if sum(set(hp_list[ins.genome_id])) == NUM_HAPLOTYPE else 'het'
            if gen_id_1[(ins.genome_id, ins.haplotype_1)]:
                db = gen_id_1[(ins.genome_id, ins.haplotype_1)][0]
                n_sup = len(set([red.read_id for red in ins.supp_read_ids]) - set(db.supp_read_ids))
                db.supp += n_sup
                db.genotype = genotype
            else:
                double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, db.bp_2, db.direction_2 ,ins.genome_id, ins.haplotype_1, ins.haplotype_1, ins.supp, ins.supp_read_ids, ins.length, genotype , 'dashed'))
                double_breaks[-1].sv_type = svtype
            
    return ins_to_remove
                    
        
def dup_to_ins(ins_list_pos, ins_list, dbs, min_sv_size, ins_clusters, double_breaks):
    
    NUM_HAPLOTYPE = 3
    db = dbs[0]
    ins_to_remove = []
    
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
        
        if db.length > ins_cl[0].length + min_sv_size:
            continue
        
        for ins in ins_cl:
            hp_list[ins.genome_id].append(ins.haplotype_1)
        
        for db in dbs:
            gen_id_1[(db.genome_id, db.haplotype_1)].append(db)
            hp_list[db.genome_id].append(db.haplotype_1)
            
        for ins in ins_cl:
            ins_to_remove.append(ins)
            genotype = 'hom' if sum(set(hp_list[ins.genome_id])) == NUM_HAPLOTYPE else 'het'
            if gen_id_1[(ins.genome_id, ins.haplotype_1)]:
                db = gen_id_1[(ins.genome_id, ins.haplotype_1)][0]
                n_sup = len(set([red.read_id for red in ins.supp_read_ids]) - set(db.supp_read_ids))
                db.supp += n_sup
                db.genotype = genotype
            else:
                double_breaks.append(DoubleBreak(db.bp_1, db.direction_1, db.bp_2, db.direction_2 ,ins.genome_id, ins.haplotype_1, ins.haplotype_1, ins.supp, ins.supp_read_ids, ins.length, genotype , 'dashed'))
                double_breaks[-1].is_dup = True
            
    return ins_to_remove
 
          
def match_long_ins(ins_clusters, double_breaks, min_sv_size):

    DEL_THR = 10000
    ins_list = defaultdict(list)
    ins_list_pos = defaultdict(list)
    ins_to_remove = []
    for ins in ins_clusters:
        ins_list_pos[ins.bp_1.ref_id].append(ins.bp_1.position)
        ins_list[ins.bp_1.ref_id].append(ins)
        
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
        
    for dbs in clusters.values():
        db = dbs[0]
        if db.bp_1.ref_id == db.bp_2.ref_id and  db.direction_1 * db.direction_2 > 0:
            continue
        if db.bp_1.ref_id == db.bp_2.ref_id and db.direction_1 > 0 and db.direction_2 < 0 and db.bp_2.position - db.bp_1.position < DEL_THR:
            continue
        if db.is_dup:
            ins_to_remove += dup_to_ins(ins_list_pos, ins_list, dbs, min_sv_size, ins_clusters, double_breaks)
        else:
            ins_to_remove_tra = tra_to_ins(ins_list_pos, ins_list, db.bp_1, db.direction_1, dbs, ins_clusters, double_breaks, min_sv_size)
            if not ins_to_remove_tra:
                ins_to_remove += tra_to_ins(ins_list_pos, ins_list, db.bp_2, db.direction_2, dbs, ins_clusters, double_breaks, min_sv_size)
            else:
                ins_to_remove += ins_to_remove_tra
                
    if ins_to_remove:
        for ins in list(set(ins_to_remove)):
            ins_clusters.remove(ins) 

def calc_vaf(db_list):
    NUM_HAPLOTYPES = 3
    for db1 in db_list.values():
        db = db1[0]
        span_bp1 = 0
        span_bp2 = 0
        for i in range(NUM_HAPLOTYPES):
            span_bp1 += db.bp_1.spanning_reads[(db.genome_id, i)]
            span_bp2 += db.bp_2.spanning_reads[(db.genome_id, i)]
        DV = 0
        DR = int(np.mean([span_bp1, span_bp2]))
        for db in db1:
            DR1 = int(np.median([db.bp_1.spanning_reads[(db.genome_id, db.haplotype_1)], db.bp_2.spanning_reads[(db.genome_id, db.haplotype_2)]]))
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
    if not control_id in sample_ids or db_list[control_id][0].vaf < control_vaf:
        mut_type = 'somatic'
    for db1 in db_list.values():
        pass_list = [db.is_pass for db in db1]
        if 'FAIL_LOWCOV_OTHER' in pass_list and not 'PASS' in pass_list:
            mut_type = 'germline'
        for db in db1:
            db.mut_type = mut_type
        
    
def annotate_mut_type(double_breaks, control_id, control_vaf):
    clusters = defaultdict(list) 
    for br in double_breaks:
        clusters[br.to_string()].append(br)
        
    for db_clust in clusters.values():
        db_list = defaultdict(list)
        for db in db_clust:
            db_list[db.genome_id].append(db)
            
        calc_vaf(db_list) 
        if control_id:
            add_mut_type(db_list, control_id, control_vaf)
        
        

def filter_germline_db(double_breaks):
    db_list = defaultdict(list) 
    for db in double_breaks:
        db_list['germline'].append(db)
        if db.mut_type == 'somatic':
            db_list['somatic'].append(db)
    return db_list

def filter_fail_double_db(double_breaks, output_only_pass, keep_low_coverage, vaf_thr):
    db_list = []
    if output_only_pass:
        for db in double_breaks:
            if db.is_pass == 'PASS' and db.vaf > vaf_thr:
                db_list.append(db)
    elif not keep_low_coverage:
        for db in double_breaks:
            if not db.is_pass == 'FAIL_LOWCOV_OTHER' and db.vaf > vaf_thr:
                db_list.append(db)
    else:
        db_list = double_breaks
    
    db_list = filter_germline_db(db_list)
    return db_list
            
    
def compute_bp_coverage(double_breaks, coverage_histograms, genome_ids):
    haplotype_list = [0, 1, 2]
    for db in double_breaks:
        if not db.bp_1.is_insertion:
            for genome_id in genome_ids:
                for haplotype in haplotype_list:
                    db.bp_1.spanning_reads[(genome_id, haplotype)] = bp_coverage(db.bp_1, genome_id, haplotype, coverage_histograms)
        if not db.bp_2.is_insertion:
            for genome_id in genome_ids:
                for haplotype in haplotype_list:
                    db.bp_2.spanning_reads[(genome_id, haplotype)] = bp_coverage(db.bp_2, genome_id, haplotype, coverage_histograms)       

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

    phased_lengths = []
    for (chr_id, block_name), coords in haplotype_blocks.items():
        if max(coords) - min(coords) > MIN_BLOCK_LEN and len(coords) >= MIN_SNP:
            endpoint_list[chr_id].append(min(coords))
            endpoint_list[chr_id].append(max(coords))
            phased_lengths.append(max(coords) - min(coords))

    for chr_id, values in endpoint_list.items():
        values.sort()
        switch_points[chr_id] = [(a + b) // 2 for a, b in zip(values[:-1], values[1:])]

    total_phased = sum(phased_lengths)
    _l50, n50 = _calc_nx(phased_lengths, total_phased, 0.50) 
    logger.info(f"\tTotal phased length: {total_phased}")
    logger.info(f"\tPhase blocks N50: {n50}")

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


def get_insertionreads(segments_by_read):
    ins_list_all = defaultdict(list)
    for read in segments_by_read.values():
        for seg in read:
            if seg.is_insertion:
                ins_list_all[seg.ref_id].append(seg)
    return ins_list_all    

def get_splitreads(segments_by_read):
    split_reads = []
    for read in segments_by_read.values():
        split = [seg for seg in read if not seg.is_insertion and not seg.is_clipped]
        if len(split)>1:
            split_reads.append(split)
    return split_reads

def resolve_overlaps(split_reads, min_ovlp_len, min_mapq):
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
            
    def _update_ovlp_seg(seglist_1, seglist_2, left_ovlp, min_mapq):
        seg_to_remove = []
        
        #if not seglist_1[-1].mapq > min_mapq and seglist_2[-1].mapq > min_mapq:
        #    seg2 = seglist_1 
        #else:
        #    seg2 = seglist_2 if seglist_2[0].strand == 1 else seglist_1
        
        seg2 = seglist_2 if seglist_2[0].strand == 1 else seglist_1
        
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

    for read_segments in split_reads:
    
        segs_to_remove =[]
        cur_cluster = []
        clusters = []
        read_segments.sort(key=lambda s:(s.align_start, s.read_start))
        
        for seg in read_segments:
            if cur_cluster and (not seg.ref_id == cur_cluster[-1].ref_id or not seg.align_start == cur_cluster[-1].align_start): 
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
                seg_to_remove = _update_ovlp_seg(clusters[i-1], clusters[i], left_ovlp, min_mapq)
                segs_to_remove += seg_to_remove
                
        for seg in list(set(segs_to_remove)):
            read_segments.remove(seg)
            
        if len(read_segments) < 2:
            split_reads.remove(read_segments)

def add_secondary_ins(double_breaks):
    for ins in double_breaks:
        if not ins.bp_2.is_insertion:
            continue
        bp_2 = Breakpoint(ins.bp_1.ref_id, ins.bp_1.position,1, ins.bp_1.qual)
        ins_2 = DoubleBreak(ins.bp_2, 1, bp_2, 1, ins.genome_id, ins.haplotype_1, ins.haplotype_2, ins.supp, ins.supp_read_ids, ins.length, ins.genotype, 'dashed')
        ins_2.is_pass = ins.is_pass
        ins_2.ins_seq = ins.ins_seq
        double_breaks.append(ins_2)
        
def cluster_db(double_breaks2):
    clusters = defaultdict(list)
    subgraph_id_list = defaultdict(list)
    for br in double_breaks2:
        br.subgraph_id = []
        if br.bp_1.is_insertion:
            continue
        clusters[br.to_string()].append(br)
    conn_inversions(clusters, subgraph_id_list)
    conn_duplications(clusters, subgraph_id_list)
    conn_inter(clusters, subgraph_id_list)
    add_subgraph_id(clusters, subgraph_id_list)

def add_subgraph_id(clusters, subgraph_id_list):
    subgraph_id = 0
    new_subgr_list = defaultdict(list)
    
    for subgr in subgraph_id_list.values():
        subgr_merge = []
        for db in subgr:
            subgr_merge += db.subgraph_id
        for id_s in list(set(subgr_merge)):
            new_subgr_list[subgraph_id] += subgraph_id_list[id_s]
        new_subgr_list[subgraph_id] = list(set(new_subgr_list[subgraph_id]))
        subgraph_id += 1
        
    pos_list = defaultdict(list)
    for subgr_id, subgr_list in new_subgr_list.items():
        pos = defaultdict(list)
        for db in subgr_list:
            pos[db.bp_1.ref_id].append(db.bp_1.position)
            pos[db.bp_2.ref_id].append(db.bp_2.position)
        for seq, pos1 in pos.items():
            pos_list[seq].append((subgr_id, min(pos1),max(pos1)))
            
    for cl in clusters.values():
        gr_flag = False
        db = cl[0]
        gr_list = pos_list[db.bp_1.ref_id]
        
        for gr in gr_list:
            if not gr_flag and db.bp_1.position > gr[1] and db.bp_1.position < gr[2]:
                gr_flag = True
                for db2 in cl:
                    db2.subgraph_id = gr[0]
                    
        if gr_flag:
            continue
        
        gr_list = pos_list[db.bp_2.ref_id]            
        for gr in gr_list:
            if not gr_flag and db.bp_2.position > gr[1] and db.bp_2.position < gr[2]:
                gr_flag = True
                for db2 in cl:
                    db2.subgraph_id = gr[0]
        
 
def conn_inversions(clusters,subgraph_id_list):
    inv_list_pos = defaultdict(list)
    inv_list_neg = defaultdict(list)
    subgraph_id = 0
    FOLD_BACK_THR = 5000
    INV_THR = 1000
    
    for cl in clusters.values():
        db = cl[0]
        if db.direction_1 == db.direction_2 == 1:
            inv_list_pos[db.bp_1.ref_id].append(cl)
            if not db.bp_2.ref_id == db.bp_1.ref_id:
                inv_list_pos[db.bp_2.ref_id].append(cl)
        elif db.direction_1 == db.direction_2 == -1:
            inv_list_neg[db.bp_1.ref_id].append(cl)
            if not db.bp_2.ref_id == db.bp_1.ref_id:
                inv_list_neg[db.bp_2.ref_id].append(cl)
                
    for seq, pos_lis in inv_list_pos.items():
        neg_lis = inv_list_neg[seq]
        pos_lis_loc = []
        neg_lis_loc = []
        pos_lis_len = []
        neg_lis_len = []
        for cl in pos_lis:
            a1 = cl[0].bp_1.position if cl[0].bp_1.ref_id == seq else cl[0].bp_2.position
            pos_lis_loc.append(a1)
            pos_lis_len.append(cl[0].bp_2.position - cl[0].bp_1.position if cl[0].bp_1.ref_id == seq else -1)
            
        for cl in neg_lis:
            a1 = cl[0].bp_1.position if cl[0].bp_1.ref_id == seq else cl[0].bp_2.position
            neg_lis_loc.append(a1)
            neg_lis_len.append(cl[0].bp_2.position - cl[0].bp_1.position if cl[0].bp_1.ref_id == seq else -1)
            
        if not pos_lis_loc or not neg_lis_loc:
            continue
        
        pos_lis_loc = np.array(pos_lis_loc)
        neg_lis_loc = np.array(neg_lis_loc)
        pr_dist = abs(pos_lis_loc[:, None] - neg_lis_loc)
        dist = pr_dist.min(axis=1)
        mx = pr_dist.max() + 1
        inv_pairs = []
        for ind, inv_len in enumerate(pos_lis_len):
            if 0 < inv_len < FOLD_BACK_THR and dist[ind] > FOLD_BACK_THR:
                for db in pos_lis[ind]:
                    db.sv_type = 'foldback'
                    pr_dist[ind] = mx
                    inv_pairs.append((ind,-1))
                    
        dist = pr_dist.min(axis=0)
        for ind, inv_len in enumerate(neg_lis_len):
            if 0 < inv_len < FOLD_BACK_THR and dist[ind] > FOLD_BACK_THR:
                for db in neg_lis[ind]:
                    db.sv_type = 'foldback'
                    pr_dist[:, ind] = mx
                    inv_pairs.append((-1,ind))
                    
        while len(inv_pairs) < min([len(pos_lis_loc) , len(neg_lis_loc)]):
            indexes = pr_dist.argmin(axis=1)
            dist = pr_dist.min(axis=1)
            pair_list = list(enumerate(indexes))
            inv_pairs.append(pair_list[0])
            for pr in pair_list[1:]:
                if not inv_pairs[-1][1] == pr[1]:
                    inv_pairs.append(pr)
                elif inv_pairs[-1][1] == pr[1] and dist[inv_pairs[-1][0]] > dist[pr[0]]:
                    inv_pairs[-1] = pr
            for pr in inv_pairs:
                pr_dist[pr[0]][pr[1]] = mx
                
        for pr in inv_pairs:
            if pr[0] == -1 or pr[0] == -1:
                continue
            clp = pos_lis[pr[0]]
            cln = neg_lis[pr[1]]
            sv_type = None
            if abs(clp[0].bp_1.position - cln[0].bp_2.position) + abs(clp[0].bp_2.position - cln[0].bp_1.position) < INV_THR:
                sv_type = 'complete_inv'
            for db in clp + cln:
                db.subgraph_id.append(subgraph_id)
                db.sv_type = sv_type
                subgraph_id_list[subgraph_id].append(db)
            subgraph_id += 1
            
def conn_duplications(clusters, subgraph_id_list):        
    
    DUP_COV_THR = 0.3
    DUP_LEN_THR = 2000000
    subgraph_id = list(subgraph_id_list.keys())[-1]+1 if len(subgraph_id_list.keys()) > 0 else 0
    for ind, cl in enumerate(clusters.values()):
        db = cl[0]
        subgraph_ls = []
        is_dup = False
        del_len = 0
        if not db.bp_1.ref_id == db.bp_2.ref_id or db.is_dup or db.bp_2.is_insertion:
            continue
        
        if db.direction_1 == -1 and db.direction_2 == 1:
            dup_len = db.bp_2.position - db.bp_1.position
            if dup_len > DUP_LEN_THR:
                continue
            
            subgraph_ls.append(cl)
            for cl2 in list(clusters.values())[ind+1:]:
                db2 = cl2[0]
                if not db2.bp_1.ref_id == db.bp_1.ref_id:
                    break
                if db2.bp_1.ref_id == db2.bp_2.ref_id and db2.bp_1.position > db.bp_1.position:
                    if not subgraph_ls:
                        is_dup = True
                    break
                if db2.bp_1.ref_id == db2.bp_2.ref_id and db.direction_1 == 1 and db.direction_2 == -1:
                    subgraph_ls.append(cl2)
                    del_len += min([db2.bp_2.position,db.bp_2.position]) - db2.bp_1.position
                else:
                    subgraph_ls.append(cl2) 
                    
            subgraph_ls.append(cl)
            subgraph_id +=1
            if not del_len:
                is_dup = True 
            elif del_len / dup_len < DUP_COV_THR:
                is_dup = True
            svtype = '' if is_dup else 'Intra_chr_ins'
            for db in cl:
                db.is_dup = is_dup
                db.sv_type = svtype
            for cl2 in subgraph_ls:
                for db2 in cl2:
                    db2.subgraph_id.append(subgraph_id)
                    subgraph_id_list[subgraph_id].append(db2)
        
            
def conn_inter(clusters, subgraph_id_list):
    subgraph_id = list(subgraph_id_list.keys())[-1]+1 if len(subgraph_id_list.keys()) > 0 else 0
    intra_chr_cl = defaultdict(list)
    DEL_THR = 1000000
    
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
                    
                    for db in cl:
                        db.subgraph_id.append(subgraph_id)
                        subgraph_id_list[subgraph_id].append(db)
                        db.sv_type = 'Inter_chr_ins'
                        
                    for db2 in cl2:
                        db2.subgraph_id.append(subgraph_id)
                        subgraph_id_list[subgraph_id].append(db2)
                        db2.sv_type = 'Inter_chr_ins'
                    subgraph_id +=1

def write_alignments(allsegments, outpath):
    aln_dump_stream = open(outpath, "w")
    for read in allsegments.values():
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
            summary_csv[br.to_string()][loc[idd]] = (br.is_pass, br.supp, br.bp_1.spanning_reads[idd], br.bp_2.spanning_reads[idd])            
        else:
            idd=(br.genome_id, br.haplotype_1)
            summary_csv[br.to_string()][loc[idd]] = (br.is_pass, br.supp, br.bp_1.spanning_reads[idd], br.bp_2.spanning_reads[idd])
            
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
            
def call_breakpoints(segments_by_read, ref_lengths, coverage_histograms, genome_ids, control_id, args):
    
    if args.write_alignments:
        outpath_alignments = os.path.join(args.out_dir, "read_alignments")
        write_alignments(segments_by_read, outpath_alignments)
        
    logger.info('Extracting split alignments')
    split_reads = get_splitreads(segments_by_read)
    logger.info('Resolving overlaps')
    resolve_overlaps(split_reads,  args.sv_size, args.min_mapping_quality)
    
    ins_list_all = get_insertionreads(segments_by_read)
    
    logger.info('Extracting clipped reads')

    clipped_clusters = []
    extract_clipped_end(segments_by_read)
    clipped_reads = get_clipped_reads(segments_by_read)
    clipped_clusters = cluster_clipped_ends(clipped_reads, args.bp_cluster_size,args.min_ref_flank, ref_lengths)
    logger.info('Starting breakpoint detection')
    double_breaks = get_breakpoints(split_reads, ref_lengths, args)
    
    logger.info('Clustering unmapped insertions')
    ins_clusters = extract_insertions(ins_list_all, clipped_clusters, ref_lengths, args)
    
    match_long_ins(ins_clusters, double_breaks, args.min_sv_size)
    
    compute_bp_coverage(double_breaks, coverage_histograms, genome_ids)
    double_breaks = double_breaks_filter(double_breaks, args.bp_min_support, genome_ids)
    double_breaks.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position, b.direction_1))
    
    compute_bp_coverage(ins_clusters, coverage_histograms, genome_ids)
    ins_clusters = insertion_filter(ins_clusters, args.bp_min_support, genome_ids)
    ins_clusters.sort(key=lambda b:(b.bp_1.ref_id, b.bp_1.position))
   
    double_breaks +=  ins_clusters
    
    if args.inbetween_ins:
        add_breakends(double_breaks, clipped_clusters, args.bp_min_support)
    
    cont_id  = list(control_id)[0] if control_id else '' 
    annotate_mut_type(double_breaks, cont_id, args.control_vaf)
        
    logger.info('Writing breakpoints')
    output_breaks(double_breaks, genome_ids, args.phase_vcf, open(os.path.join(args.out_dir,"breakpoints_double.csv"), "w"))
    
    double_breaks = filter_fail_double_db(double_breaks, args.output_only_pass, args.keep_low_coverage, args.vaf_thr)
    
    return double_breaks 
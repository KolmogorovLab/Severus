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
    __slots__ = "ref_id", "position", "spanning_reads", "connections" , "read_ids", "pos2"
    def __init__(self, ref_id, ref_position):
        self.ref_id = ref_id
        self.position = ref_position
        self.spanning_reads = defaultdict(int)
        self.connections =[]
        self.read_ids=[]
        self.pos2 = []


class DoubleBreak(object):
    __slots__ = ("bp_1", "direction_1", "bp_2", "direction_2", "genome_id","haplotype1",'haplotype2',
                 "supp",'supp_read_ids','length','genotype','edgestyle','bp_id')
    def __init__(self, bp_1, direction_1, bp_2, direction_2, genome_id,haplotype1,haplotype2,
                 supp,supp_read_ids,length, genotype, edgestyle,bp_id):
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
        self.bp_id = bp_id
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
    __slots__ = "genome_id","haplotype", "ref_id", "dir1", "pos1", "dir2" , "pos2", "coverage","length_bp"
    def __init__(self, genome_id,haplotype,ref_id, pos1, pos2,coverage,length_bp):
        self.genome_id=genome_id
        self.haplotype=haplotype
        self.ref_id = ref_id
        self.dir1 = '-'
        self.pos1 = pos1
        self.dir2 = '+'
        self.pos2 = pos2
        self.coverage = coverage
        self.length_bp=length_bp

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

####AYSE: CHECK VCF FROM OTHER TOOLS
def get_phasingblocks(hb_vcf):
    vccf=pysam.VariantFile(hb_vcf)
    hb = defaultdict(list)
    hb_list = defaultdict(list)
    hb_points = defaultdict(list)
    for var in vccf:
        if 'PS' in var.samples.items()[0][1].items()[-1]:
            hb[(var.chrom, var.samples.items()[0][1]['PS'])].append(var.pos)
    for key,hb1 in hb.items():
        hb_list[key[0]].append(min(hb1))
        hb_list[key[0]].append(max(hb1))
    for key, values in hb_list.items():
        values.sort()
        hb_points[key]=([int((a + b)/2) for a, b in zip(values[1::2], values[2::2])])
    return hb_points

def get_genomicsegments(double_breaks,bam_files, thread_pool,hb_vcf):
    if hb_vcf:
        hb_points=get_phasingblocks(hb_vcf)
        phased=True
    else:
        hb_points = defaultdict(list)
        phased = False
    single_bp = defaultdict(list)
    genomicsegments=[]
    segments=[]
    for double_bp in double_breaks:        
        single_bp[(double_bp.genome_id ,double_bp.haplotype1, double_bp.bp_1.ref_id)].append(double_bp.bp_1.position)
        single_bp[(double_bp.genome_id ,double_bp.haplotype2, double_bp.bp_2.ref_id)].append(double_bp.bp_2.position)
    for key,s_bp in single_bp.items():
        s_bp1= s_bp+hb_points[key[2]]
        s_bp1 = list(set(s_bp1))
        s_bp1.sort()
        for a, b in zip(s_bp1[:-1], s_bp1[1:]):
            bam = [bbam for bbam in bam_files if key[0] in bbam]
            segments.append((bam,key[2], a, b,key[1]))
    genomicsegments = get_segments_coverage(segments, thread_pool,phased)
    return genomicsegments , hb_points
                        
  

##Ayse: haplotype info/ single_bp
def get_breakpoints(allreads_bisect, split_reads, clust_len, max_unaligned_len,
                    min_reads, min_ref_flank, ref_lengths, min_mapq,single_bp,lowmapq_reg,min_sv_size,ins_list_new):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """
    seq_breakpoints = defaultdict(list)
    print('breakpoint_finder_1')

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
        seq_breakpoints[s2.ref_id].append(ReadConnection(s2.ref_id, ref_bp_2, sign_2, s1.ref_id, ref_bp_1, sign_1,
                                                             s2.haplotype, s1.haplotype, s2.read_id, s2.genome_id))
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
            pos2 = []
            for x in cl:
                unique_reads.add((x.read_id, (x.genome_id,x.haplotype_1)))
                read_ids.append(x.read_id)
                pos2.append(x.signed_coord_2())
            by_genome_id = defaultdict(int)
            for read in unique_reads:
                by_genome_id[read[1]] += 1
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.pos_1 for x in cl]))
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and lowmapq and lowmapq[1][bisect.bisect_left(lowmapq[0],position)-1]<position:
                    bp_cluster = Breakpoint(seq, position)
                    bp_cluster.connections = cl
                    bp_cluster.read_ids = read_ids
                    bp_cluster.pos2=pos2
                    bp_clusters[seq].append(bp_cluster)
                    
    for key, bp_cluster in bp_clusters.items():####sort with end as well!!
        read_segments = allreads_bisect[key]
        read_segments[0], read_segments[1], read_segments[2], read_segments[3] = zip(*sorted(zip(read_segments[0],read_segments[1],read_segments[2],read_segments[3])))
        for bp in bp_cluster:
            strt=bisect.bisect_left(read_segments[1],bp.position)
            end = bisect.bisect_left(read_segments[0],bp.position)
            count_all = Counter([read_segments[3][i] for i in range(strt,end) if read_segments[1][i] > bp.position])
            for gen_id,counts in count_all.items():
                bp.spanning_reads[gen_id]=counts
    double_breaks = get_2_breaks(bp_clusters, clust_len, min_reads,min_sv_size)
    ins_clusters = extract_insertions(ins_list_new, lowmapq_reg,clust_len,min_ref_flank,ref_lengths,min_reads,allreads_bisect)
    double_breaks +=  ins_clusters             
    double_breaks.sort(key=lambda b: (b.bp_1.ref_id, b.bp_1.position, b.direction_1))

    
    return double_breaks

## Fixed the _normalize_coord part with new clustering 
def get_2_breaks(bp_clusters, clust_len, min_reads,min_sv_size):
    def _normalize_coord(read_id, clusters, coord):
        for cl1 in clusters:
            if read_id in cl1.read_ids and coord in cl1.pos2:
                return int(math.copysign(1, coord)), cl1
        return None, None
    double_breaks = []
    bp_num = 1
    double_connections = defaultdict(list)
    chlist = list(bp_clusters.keys())
    for bp_cluster in bp_clusters.values():
        for cl in bp_cluster:
            if len(cl.connections)<1000:## Ask Misha
                for conn in cl.connections:
                    if conn.ref_id_2 is None:
                        continue
                    if conn.ref_id_2 not in chlist:
                        continue
                    dir_2, bp_2 = _normalize_coord(conn.read_id, bp_clusters[conn.ref_id_2], conn.signed_coord_1())
                    dir_1, bp_1 = _normalize_coord(conn.read_id, bp_clusters[conn.ref_id_1], conn.signed_coord_2())
                    if None in [bp_1, bp_2]:
                        continue
                    if bp_2.position - bp_1.position > min_sv_size:
                        double_connections[(bp_1, dir_1, bp_2, dir_2)].append(conn)
    for (bp_1, dir_1, bp_2, dir_2), conn_list in double_connections.items():
        support = defaultdict(int)
        support_reads = defaultdict(list)
        hap1_support= defaultdict(int)
        hap2_support= defaultdict(int)
        for conn1 in conn_list:
            support[(conn1.genome_id, conn1.haplotype_1,conn1.haplotype_2)] +=1
            support_reads[(conn1.genome_id, conn1.haplotype_1,conn1.haplotype_2)]+=[conn1.read_id]
            hap1_support[conn1.haplotype_1]+=1
            hap2_support[conn1.haplotype_1]+=1            
        if max(support.values())>= min_reads:
            for (genome_id,haplotype1,haplotype2),supp in support.items():
                if supp>0:
                    if (hap1_support[1] and hap1_support[2]) or (hap2_support[1] and hap2_support[2]):
                        genotype = "hom"
                    else:
                        genotype = "het"
                    if bp_1.ref_id == bp_1.ref_id:
                        length_bp = abs(bp_1.position - bp_2.position)
                    else:
                        length_bp = 0
                    bp_id = 'bp_'+str(bp_num)
                    bp_num+=1
                    double_breaks.append(DoubleBreak(bp_1, dir_1, bp_2, dir_2,genome_id,haplotype1, haplotype2, supp,
                                                     support_reads[(genome_id,haplotype1,haplotype2)], length_bp, genotype, 'dashed', bp_id))
    
    return double_breaks


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
        bp = ','.join(bp_array)
        bp_to_write = ','.join([key, ','.join(bp_array)])
        out_stream.write(bp_to_write)
        out_stream.write("\n")
        
def extract_insertions(ins_list_new, lowmapq_reg,clust_len,min_ref_flank,ref_lengths,min_reads,allreads_bisect):
    ins_list = defaultdict(list)
    for ins in ins_list_new:
        ins_list[ins[0]].append(ins)
    ins_clusters = []
    t = 0
    for seq, ins_pos in ins_list.items():
        clusters = []
        cur_cluster = []
        ins_pos.sort(key=lambda x:(x[1]))
        lowmapq = lowmapq_reg[seq]
        read_segments = allreads_bisect[seq]
        read_segments[0], read_segments[1], read_segments[2], read_segments[3] = zip(*sorted(zip(read_segments[0],read_segments[1],read_segments[2],read_segments[3])))
        for rc in ins_pos:
            if cur_cluster and rc[1] - cur_cluster[-1][1] > clust_len:
                cur_cluster.sort(key=lambda x:x[5])
                cl_ins = []
                for cl1 in cur_cluster:
                    if cl_ins and abs(cl1[5] - cl_ins[-1][5]) > clust_len:
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
            unique_reads = set()
            read_ids = []
            pos2 = []
            for x in cl:
                unique_reads.add((x[2], (x[3],x[4])))
                read_ids.append(x[2])
                pos2.append(x[5])
            by_genome_id = defaultdict(list)
            for read in unique_reads:
                by_genome_id[read[1]].append(read[0])
            if max([len(value) for value in by_genome_id.values()]) >= min_reads:
                position = int(np.median([x[1] for x in cl]))
                length = int(np.median([x[5] for x in cl]))
                edgeW = defaultdict(int)
                for key in by_genome_id.keys():
                    edgeW[key[0]]+=1
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank and lowmapq and lowmapq[1][bisect.bisect_left(lowmapq[0],position)-1]<position:
                    ins_id = 'INS_'+str(t)
                    t+=1
                    for key, value in by_genome_id.items():
                        bp_1 = Breakpoint(seq, position)
                        bp_2 = Breakpoint(seq, position)
                        bp_1.read_ids = value
                        bp_2.read_ids = value
                        count_all = Counter(read_segments[3][bisect.bisect_left(read_segments[1],position):bisect.bisect_left(read_segments[0],position)])
                        for gen_id,counts in count_all.items():
                            bp_1.spanning_reads[gen_id]=counts
                            bp_2.spanning_reads[gen_id]=counts

                        genotype = "hom" if len(edgeW) > 1 else "het"
                        db = DoubleBreak(bp_1, -1, bp_2, 1, key[0], key[1], key[1], len(value), value, length, genotype, 'dashed', ins_id)
                        ins_clusters.append(db)
    return(ins_clusters)
                    
     
        
def get_segments_coverage(segments, thread_pool,phased):
    print("Computing segments coverage", file=sys.stderr)
    random.shuffle(segments)
    tasks = [(s[0][0], s[1], s[2],s[3],s[4],phased) for s in segments]
    seg_coverage = None
    seg_coverage = thread_pool.starmap(_get_median_depth, tasks)
    genomicsegments = []
    for seg, coverage in zip(segments, seg_coverage):
        genome_id = os.path.basename(seg[0][0])
        genomicsegments.append(GenomicSegment(genome_id, seg[4], seg[1], seg[2], seg[3],coverage,abs(seg[3]-seg[2])))
    return genomicsegments

def _get_median_depth(bam, ref_id, ref_start, ref_end,hp,phased):
    SAMTOOLS_BIN = "/data/KolmogorovLab/tahmad/tools/samtools/samtools"
    cov_by_bam = []
    try:
        if not phased:
            samtools_out = subprocess.Popen("{0} coverage {1} -r '{2}:{3}-{4}' -q 10 -l 100 -d 1000"
                                            .format(SAMTOOLS_BIN, bam, ref_id, ref_start, ref_end),
                                            shell=True, stdout=subprocess.PIPE).stdout
        else:
            samtools_out = subprocess.Popen("{0} coverage {1} -r '{2}:{3}-{4}' -p {5} -q 10 -l 100 -d 1000"
                                            .format(SAMTOOLS_BIN, bam, ref_id, ref_start, ref_end, hp),
                                            shell=True, stdout=subprocess.PIPE).stdout
        for line in samtools_out:
            if line.startswith(b"#"):
                continue
            fields = line.split()
            coverage = float(fields[6])
            cov_by_bam.append(coverage)
            break
    except Exception as e:
        print(e)
        print("{0} coverage {1} -r '{2}:{3}-{4}' -p {5} -q 10 -l 100 -d 1000"
              .format(SAMTOOLS_BIN, bam, ref_id, ref_start, ref_end , hp))
    if len(cov_by_bam) > 0:
        return cov_by_bam
    else:
        return 0
        


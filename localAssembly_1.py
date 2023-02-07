#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:11:59 2023

@author: keskusa2
"""
import subprocess
import pysam
import os
import shutil
import numpy as np
from collections import  defaultdict
from multiprocessing import Pool
from multiprocessing import Pool
from breakpoint_finder import resolve_overlaps
from bam_processing import ReadSegment

class db_cluster(object):
    __slots__ = "genome_id","haplotype", "pos", "read_ids", 'bp_list'
    def __init__(self, genome_id,haplotype,pos, read_ids,bp_list):
        self.genome_id=genome_id
        self.haplotype=haplotype
        self.pos = pos
        self.read_ids = read_ids
        self.bp_list=bp_list
        
class Breakpoint(object):
    __slots__ ="read_start", "ref_id", "position", "spanning_reads"
    def __init__(self, ref_id, ref_position):
        self.ref_id = ref_id
        self.position = ref_position
        self.spanning_reads = defaultdict(int)


class DoubleBreak(object):
    __slots__ = "bp_1", "direction_1", "bp_2", "direction_2", "genome_id","haplotype1",'haplotype2',"supp",'length','edgeweight'
    def __init__(self, bp_1, direction_1, bp_2, direction_2, genome_id,haplotype1,haplotype2,supp,length,edgeweight):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
        self.direction_2 = direction_2
        self.genome_id = genome_id
        self.haplotype1 = haplotype1
        self.haplotype2 = haplotype2
        self.supp = supp
        self.length = length
        self.edgeweight = edgeweight
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


def extract_reads_from_bam(all_bams, outpath,db, asm_id):
    out_folder = outpath+asm_id
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    bam_file = [bbam for bbam in all_bams if db.genome_id in bbam]
    aln_file = pysam.AlignmentFile(bam_file[0], "rb")
    out_file = out_folder+'/extracted.fa'
    written_reads = []
    fast=open(out_file,'w')
    for pos1 in db.pos:
        aln = aln_file.fetch(pos1[0],pos1[1]-100, pos1[1]+100)
        for read in aln:
            if read.qname in db.read_ids and not read.qname in written_reads and not read.is_supplementary and not read.is_secondary:
                fast.write('>'+read.qname+'\n'+read.seq+'\n')
                written_reads.append(read.qname)
    fast.close()
        
            
def call_flye(outpath, platform, nthreads):
    fasta_file = outpath+'/extracted.fa'
    if platform == 'ont':
        flye_arg = '--nano-raw'
    elif platform == 'pacbio':
        flye_arg = '--pacbio-hifi'
    subprocess.run("flye {0} {1} --out-dir {2} --threads {3} --min-overlap 1000 --meta --genome-size 30k"
                   .format(flye_arg,fasta_file,outpath,nthreads),shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
def call_minimap2(flye_out, ref_fasta,nthreads,out_bam):
    subprocess.run("minimap2 -ax asm5 -k 17 -K 5G {0} {1} -t {2} | samtools sort -@{2} -m 4G > {3}"
                   .format(ref_fasta,flye_out, nthreads,out_bam),shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    subprocess.run("samtools index -@ 4 {0}".format(out_bam), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
def localAssembly(db_id, outpath_assembly,platform,nthreads,ref_fasta):
    outpath = outpath_assembly+db_id
    bam_list= []
    flye_out = outpath+'/assembly.fasta'
    out_bam = outpath_assembly+'/bam/'+db_id+'.bam'
    call_flye(outpath, platform, nthreads)
    if os.path.exists(flye_out):
        call_minimap2(flye_out, ref_fasta, nthreads,out_bam)
        bam_list.append(os.path.basename(flye_out)[:-4])
        shutil.rmtree(outpath)
    return bam_list

##combine with bamprocessing 

def breakpoints_from_assembly(read_segments,db_list,min_svsize):
    def _signed_breakpoint(seg, direction):
        ref_bp, sign = None, None
        if direction == "right":
            ref_bp = seg.ref_end if seg.strand == "+" else seg.ref_start
            sign = 1 if seg.strand == "+" else -1
        elif direction == "left":
            ref_bp = seg.ref_start if seg.strand == "+" else seg.ref_end
            sign = -1 if seg.strand == "+" else 1
        return ref_bp, sign
    double_break =[]
    for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
        bp_1_position, direction_1 = _signed_breakpoint(s1, "right")
        bp_2_position, direction_2 = _signed_breakpoint(s2, "left")
        bp_1 = Breakpoint(s1.ref_id, bp_1_position)
        bp_2 = Breakpoint(s2.ref_id, bp_2_position)
        length=0
        dist = 100000000
        if s1.ref_id==s2.ref_id:
            length = abs(bp_1_position-bp_2_position)
        if length == 0 or length > min_svsize:
            if len(db_list)>1:
                for db1 in db_list:
                    if db1.bp_1.ref_id == s1.ref_id and db1.bp_2.ref_id == s2.ref_id:
                        dist2 = np.median([abs(bp_1_position-db1.bp_1.position),abs(bp_2_position-db1.bp_2.position)])
                        if dist2 < dist:
                            dist = dist2
                            db = db1
            if dist < 100000000:
                bp_1.spanning_reads = db.bp_1.spanning_reads
                bp_2.spanning_reads = db.bp_2.spanning_reads
                double_break.append(DoubleBreak(bp_1, direction_1, bp_2, direction_2, db.genome_id, db.haplotype1, db.haplotype2, db.supp, length, db.edgeweight))
            else:
                if s1.ref_id != s2.ref_id:
                    continue
                else:
                    for db1 in db_list:
                        if db1.bp_1.ref_id == s1.ref_id:
                            bp_1.spanning_reads = db1.bp_1.spanning_reads
                            bp_2.spanning_reads = db1.bp_1.spanning_reads
                            double_break.append(DoubleBreak(bp_1, direction_1, bp_2, direction_2, db1.genome_id, db1.haplotype1, db1.haplotype1, db1.supp, length, db1.edgeweight))
                        elif  db1.bp_2.ref_id == s1.ref_id:
                            bp_1.spanning_reads = db1.bp_2.spanning_reads
                            bp_2.spanning_reads = db1.bp_2.spanning_reads
                            double_break.append(DoubleBreak(bp_1, direction_1, bp_2, direction_2, db1.genome_id, db1.haplotype2, db1.haplotype2, db1.supp, length, db1.edgeweight))                       
    return double_break
                                

def add_insertion(strand, ins_start,ins_len,ins_id,ref_id,db_list,double_break_ins):
    dist = 100000000
    for db1 in db_list:
        if db1.bp_1.ref_id == ref_id:
            dist2 = abs(ins_start-db1.bp_1.position)
            if dist2 <dist:
                dist = dist2
                db = db1
                haplotype= db.haplotype1
        elif  db1.bp_1.ref_id == ref_id:
            dist2 = abs(ins_start-db1.bp_2.position)
            if dist2 <dist:
                dist = dist2
                db = db1
                haplotype= db.haplotype2
    bp_1 = Breakpoint(ref_id, ins_start)
    bp_2=Breakpoint(ins_id, ins_len)
    double_break_ins.append(DoubleBreak(bp_1, -1, bp_2, 1, db.genome_id, haplotype, haplotype, db.supp, ins_len, db.edgeweight))
    double_break_ins.append(DoubleBreak(bp_2, -1, bp_1, -1, db.genome_id, haplotype,haplotype, db.supp, ins_len, db.edgeweight))
    return double_break_ins
    
def find_segments(outpath_bam,asm_ids, db_list, min_svsize,min_ovlp_len,sv_clustersize):
    double_break_list = defaultdict(list)
    for asm_id in asm_ids:
        bam_path = outpath_bam + asm_id+'.bam'
        aln_file = pysam.AlignmentFile(bam_path, "rb")
        db = db_list[asm_id]
        ins_list = defaultdict(str)
        read_segments= []
        double_break_ins = []
        haplotype=db[0].haplotype1
        genome_id= db[0].genome_id
        for read in aln_file:
            cigar = read.cigar
            read_start=0
            read_aligned = 0
            ref_aligned = 0 
            ref_start = read.reference_start
            strand = '+'
            read_length = np.sum([k[1] for k in cigar if k[0]!=2])
            first_clip =True
            for k in cigar:
                if k[0]>2:
                    if first_clip:
                        read_start = k[1]
                first_clip = False
                if k[0]==0: 
                    read_aligned +=k[1]
                    ref_aligned+=k[1]
                if k[0] == 1:
                    if k[1]>min_svsize:
                        ins_start = ref_start + ref_aligned
                        ins_start_read = read_start + read_aligned-1
                        ins_id = 'INS_'+read.reference_name+'_'+str(ins_start)
                        ins_list[ins_id] = read.seq[ins_start_read:ins_start_read+k[1]]
                        double_break_ins = add_insertion(strand, ins_start,k[1],ins_id,read.reference_name,db,double_break_ins)
                        read_aligned+=k[1]
                    else:
                        read_aligned+=k[1]
                if k[0] == 2:
                    if k[1]<min_svsize:
                        ref_aligned+=k[1]
                    elif k[1]>min_svsize:
                        ref_end = ref_start + ref_aligned
                        read_end = read_start + read_aligned
                        if read.is_reverse:
                            strand = '-'
                            read_start1, read_end1 = read_length - read_end, read_length - read_start
                        read_segments.append(ReadSegment(read_start1, read_end1, ref_start, ref_end, asm_id,
                                                read.reference_name, strand, read_length, haplotype, read.mapq, genome_id))
                        read_start = read_end+1
                        ref_start = ref_end+k[1]+1
                        read_aligned = 0
                        ref_aligned = 0
            ref_end = ref_start + ref_aligned
            read_end = read_start + read_aligned
            if read.is_reverse:
                strand = '-'
                read_start, read_end = read_length - read_end, read_length - read_start
            read_segments.append(ReadSegment(read_start, read_end, ref_start, ref_end, asm_id,
                                    read.reference_name, strand, read_length, haplotype, read.mapq, genome_id))
        read_segments.sort(key=lambda s: s.read_start)
        read_segments = resolve_overlaps([read_segments], min_ovlp_len)
        double_break = breakpoints_from_assembly(read_segments[0],db, min_svsize)
        double_break_list[asm_id] = double_break + double_break_ins
    first_item = True
    double_breaks=[]
    for value in double_break_list.values():
        if first_item:
            for db in value:
                double_breaks.append(db)
            first_item = False
        else:
            for db in value:
                for db1 in double_breaks:
                    if db.bp_1.ref_id == db1.bp_1.ref_id and db.bp_2.ref_id == db1.bp_2.ref_id:
                        if abs(db.bp_1.position - db1.bp_1.position) <= sv_clustersize and abs(db.bp_2.position - db1.bp_2.position) <= sv_clustersize:
                            db.bp_1.position = db1.bp_1.position
                            db.bp_2.position = db1.bp_2.position
                            break
                    elif db.bp_1.ref_id == db1.bp_1.ref_id and 'INS' in db.bp_2.ref_id:
                        if abs(db.bp_1.position - db1.bp_1.position) <= sv_clustersize and abs(db.bp_2.position - db1.bp_2.position) <= sv_clustersize:
                            db.bp_1.position = db1.bp_1.position
                            db.bp_2.position = db1.bp_2.position
                            db.bp_2.ref_id = db1.bp_2.ref_id
                            break
                    elif db.bp_2.ref_id == db1.bp_2.ref_id and 'INS' in db.bp_1.ref_id:
                        if abs(db.bp_1.position - db1.bp_1.position) <= sv_clustersize and abs(db.bp_2.position - db1.bp_2.position) <= sv_clustersize:
                            db.bp_1.position = db1.bp_1.position
                            db.bp_2.position = db1.bp_2.position
                            db.bp_1.ref_id = db1.bp_1.ref_id
                            break
                double_breaks.append(db)
    return double_breaks , ins_list
        
def localAssembly_main(all_bams, double_breaks,outpath,nthreads,platform,ref_fasta,min_reads, min_svsize, min_ovlp_len,sv_clustersize):
    dbcluster_list = defaultdict(list)
    db_list = defaultdict(list)
    outpath_assembly = outpath+'/assembly_files/'
    outpath_bam = outpath+'/assembly_files/bam/'
    if not os.path.exists(outpath_bam):
        os.makedirs(outpath_bam)
    t = 0
    for db in double_breaks:
        add_ext = []
        for key, value in dbcluster_list.items():
            add_ext = [key for read_id in db.supp_read_ids if read_id in value.read_ids]
            if add_ext:
                asm_id = add_ext[0]
                dbcluster_list[asm_id].read_ids = list(set(value.read_ids + db.supp_read_ids))
                dbcluster_list[asm_id].bp_list += [db.bp_id]
                dbcluster_list[asm_id].pos+=[(db.bp_1.ref_id ,db.bp_1.position)]
                dbcluster_list[asm_id].pos+=[(db.bp_2.ref_id ,db.bp_2.position)]
                db_list[asm_id].append(db)
                break
        if not add_ext:
            asm_id= 'asm_'+str(t)
            t+=1
            db_list[asm_id] =[db]
            dbcluster_list[asm_id] = db_cluster(db.genome_id,db.haplotype1, [(db.bp_1.ref_id ,db.bp_1.position),(db.bp_2.ref_id ,db.bp_2.position)], db.supp_read_ids,[db.bp_id])
    tasks = [(all_bams,outpath_assembly, db, asm_id) for asm_id,db in dbcluster_list.items()]
    thread_pool = Pool(nthreads)
    res=thread_pool.starmap(extract_reads_from_bam, tasks)
    nthreads_flye = 4 #### Ask Misha
    tasks_localAssembly=[(asm_id, outpath_assembly,platform,nthreads_flye,ref_fasta) for asm_id in dbcluster_list.keys()]
    nt = nthreads//nthreads_flye
    thread_pool.close()
    thread_pool.join()
    thread_pool_nt = Pool(nt) ## Ask Misha 4 or 1 thread per task
    print('Starting Local Assembly')
    res = thread_pool_nt.starmap(localAssembly, tasks_localAssembly)
    thread_pool_nt.close()
    thread_pool_nt.join()
    print('Local Assembly Done')
    thread_pool = Pool(nthreads)
    bamlist = os.listdir(outpath_bam)
    bam_list = []
    for bam in bamlist:
        if bam.endswith('.bam'):
            bam_list.append(bam[:-4])
    bamlist = cluster_assemblies(dbcluster_list , bam_list)
    print('Breakpoint detection')
    tasks = [(outpath_bam,asm_ids, db_list, min_svsize,min_ovlp_len,sv_clustersize) for asm_ids in bamlist]
    parsing_results = thread_pool.starmap(find_segments, tasks)
    double_breaks=[]
    ins_list = []
    for double_break1,ins_list1 in parsing_results:
        double_breaks += double_break1
        ins_list = merge(ins_list, ins_list1) ## write to fasta
    return double_breaks , ins_list
        
    
def cluster_assemblies(dbcluster_list , bamlist):
    curr_db=[]
    unique_db_list = []
    for asm_id,db_clust in dbcluster_list.items():
        if asm_id in bamlist:
            if not curr_db:
                curr_db = db_clust
                curr_list = [asm_id]
            elif np.sum([1 for pos in db_clust.pos if pos in curr_db.pos])>1:
                curr_list.append(asm_id)
            else:
                unique_db_list.append(curr_list)
                curr_list = [asm_id]
                curr_db = db_clust
    if curr_db:
        unique_db_list.append(curr_list)
    return unique_db_list


    
    
    
### TODO
# cases with less support than min_reads in other samples



##### Insertions
# add sofclipped reads as insertion


#calc the mismatch ratio for each segment
#compare with NM in full aligns     



###Check COLO graph















    
        
        
        
    
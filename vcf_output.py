#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pysam
from collections import namedtuple, defaultdict, Counter
from datetime import datetime
import sys
import os
import math

class vcf_format(object):
    __slots__ = ('Chrom', 'pos', 'haplotype', 'ID', 'sv_type', 'sv_len', 'qual', 'Filter', 'chr2', 'pos2', 'DR', 'DV', 'mut_type', 'hVaf')
    def __init__(self, Chrom, pos, haplotype, ID, sv_type, sv_len, qual, Filter, chr2, pos2, DR, DV, mut_type, hVaf):
        self.Chrom = Chrom
        self.pos = pos
        self.haplotype = haplotype
        self.ID = ID
        self.sv_type = sv_type 
        self.sv_len = sv_len
        self.qual = qual
        self.Filter = Filter
        self.chr2 = chr2
        self.pos2 = pos2
        self.DR = DR
        self.DV = DV
        self.hVaf = hVaf
        self.mut_type = mut_type
    def vaf(self):
        return self.DV / (self.DV + self.DR) if self.DV > 0 else 0
    def hVAF(self):
        return '{0:.2f}|{1:.2f}|{2:.2f}'.format(self.hVaf[0], self.hVaf[1], self.hVaf[2])
    def call_genotype(self):
        err = 0.1/2
        prior = 1/3
        genotypes = ["0/0", "0/1", "1/1"]
        MAX_READ = 100
        if self.DV + self.DR > MAX_READ:
            DV = int(MAX_READ*self.DV/(self.DV+self.DR))
            DR = MAX_READ - DV
        else:
            DV, DR = self.DV, self.DR
        hom_1 = float(pow((1-err), DV) * pow(err, DR) * prior)
        hom_0 = float(pow(err, DV) * pow((1-err), DR) * prior)
        het = float(pow(0.5, DV + DR) * prior)
        g_list = [hom_0, het, hom_1]
        g_list = [-10 * math.log10(p) for p in g_list]
        min_pl = min(g_list)
        g_list = [pl - min_pl for pl in g_list]
        GT = genotypes[g_list.index(max(g_list))]
        g_list.sort()
        GQ = g_list[1]
        PL = g_list[0]
        return GT, GQ, PL
    def info(self):
        return "SVLEN={0};SVTYPE={1};CHR2={2};END={3};MAPQ={4};SUPPREAD={5};hVAF={6}".format(
            self.sv_len, self.sv_type, self.chr2, self.pos2, self.qual, self.DV, self.hVAF())
    def sample(self):
        GT, GQ, PL = self.call_genotype()
        return "{0}:{1:.2f}:{2:.2f}:{3:.2f}{4}{5}".format(GT, GQ, PL, self.vaf(), self.DR, self.DV)
    def to_vcf(self):
        return "{0}\t{1}\t{2}\t N \t <{3}>\t{4}\t{5}\t{6}\t GT:GQ:PL:VAF:DR:DV\t{7}\n".format(
            self.Chrom, self.pos, self.ID, self.sv_type, self.qual, self.Filter, self.info(), self.sample())

  
 
def get_sv_type(db):
    if db.bp_2.is_insertion or db.bp_1.is_insertion:
        return 'INS'
    if not db.bp_1.ref_id == db.bp_2.ref_id:
        return 'BND'
    if db.direction_1 == db.direction_2:
        return 'INV'
    if db.bp_1.dir_1 < 0:
        return 'BND'
    if db.bp_1.dir_1 > 0:
        return 'DEL'

def db_2_vcf(double_breaks, control_id):
    NUM_HAPLOTYPES = 3
    t = 0
    vcf_list = defaultdict(list)
    clusters = defaultdict(dict) 
    for br in double_breaks:
    	if not br.bp_1.is_insertion:
            clusters[br.to_string()][br.genome_id] = br
    for db_clust in clusters.values():
        mut_type = 'germline'#
        sample_ids = list(db_clust.keys())
        if not control_id in sample_ids:
            mut_type = 'somatic'
        sv_type = get_sv_type(db_clust[sample_ids[0]])
        if not sv_type:
            continue
        db_list = defaultdict(list)
        for db in db_clust.values():
            db_list[db.genome_id].append(db)
        for db1 in db_list.values():
            hVaf = defaultdict(list)
            for i in range(NUM_HAPLOTYPES):
                hVaf[i]=0.0
            ID = 'SEVERUS_' + sv_type + str(t)
            t += 1
            DR = 0
            DV = 0
            qual_list = []
            pass_list = [db.is_pass for db in db1]
            if 'PASS' in pass_list:
                sv_pass = 'PASS' 
            elif 'PASS_LOWCOVERAGE' in pass_list:
                sv_pass = 'PASS_LOWCOVERAGE' 
                mut_type = 'germline'
            else:
                sv_pass = db1[0].is_pass
            for db in db1:
                DR1 = int(np.median([db.bp_1.spanning_reads[(db.genome_id, db.haplotype_1)], db.bp_2.spanning_reads[(db.genome_id, db.haplotype_2)]]))
                DV1 = db.supp
                hVaf[db.haplotype_1] =  DV1 / (DV1 + DR1) if DV1 > 0 else 0
                DR += DR1
                DV += DV1
                if db.is_pass == 'PASS':
                    qual_list += [db.bp_1.qual, db.bp_2.qual]
            if not qual_list:
                qual_list = [0]
            vcf_list[db.genome_id].append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, db.length, int(np.median(qual_list)), 
                                                     sv_pass, db.bp_2.ref_id, db.bp_2.position, DR, DV, mut_type, hVaf))#
    return vcf_list
            
def write_vcf_header(ref_lengths, outfile):
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write('##source=SEVERUS\n')
    outfile.write('##CommandLine= '+ " ".join(sys.argv[1:]) +'\n')
    filedate = str(datetime.now()).split(' ')[0]
    outfile.write('##fileDate='+filedate+'\n')#
    for chr_id, chr_len in ref_lengths.items():
        outfile.write("##contig=<ID={0},length={1}>\n".format(chr_id, chr_len))#
    outfile.write('##ALT=<ID=DEL,Description="Deletion">\n')
    outfile.write('##ALT=<ID=INS,Description="Insertion">\n')
    outfile.write('##ALT=<ID=DUP,Description="Duplication">\n')
    outfile.write('##ALT=<ID=INV,Description="Inversion">\n')
    outfile.write('##ALT=<ID=BND,Description="Breakend">\n')#
    outfile.write('##FILTER=<ID=PASS,Description="All filters passed">\n') ## add other filters as well!!#
    outfile.write("##INFO=<ID=HAPLOTYPE,Number=1,Type=String,Description=\"Haplotype of the SV\">\n")
    outfile.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    outfile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    outfile.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">\n")
    outfile.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV\">\n")
    outfile.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of supporting reads\">\n")
    outfile.write("##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description=\"Number of supporting reads\">\n")#
    outfile.write("##INFO=<ID=hVAF,Number=1,Type=String,Description=\"# Haplotype specific variant Allele frequency (H0|H1|H2)\">\n")
    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    outfile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotyping Quality\">\n")
    outfile.write("##FORMAT=<ID=PL,Number=1,Type=Integer,Description=\"PL\">\n")
    outfile.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# Number of spanning reads\">\n")
    outfile.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# Number of supporting reads\">\n")
    outfile.write("##FORMAT=<ID=VAF,Number=1,Type=Integer,Description=\"# Variant Allele frequency\">\n")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    

def write_somatic_vcf(vcf_list, outfile):
    for db in vcf_list:
        if db.mut_type == 'somatic':
            outfile.write(db.to_vcf())
    outfile.close()
    
def write_germline_vcf(vcf_list, outfile):
    for db in vcf_list:
        outfile.write(db.to_vcf())
    outfile.close()
    
    
def write_to_vcf(double_breaks, target_ids, control_id, outpath, ref_lengths):
    control_id = list(control_id)
    if not control_id:
        control_id = ['']
    vcf_list = db_2_vcf(double_breaks, control_id[0])
    if control_id:
        for target_id in target_ids:
            somatic_outfile = open(os.path.join(outpath,"SEVERUS_somatic_" + target_id + ".vcf"), "w")
            write_vcf_header(ref_lengths,somatic_outfile)
            write_somatic_vcf(vcf_list[target_id], somatic_outfile)
    all_ids = list(target_ids) + control_id
    for target_id in all_ids:
        germline_outfile = open(os.path.join(outpath,"SEVERUS_" + target_id + ".vcf"), "w")
        write_vcf_header(ref_lengths, germline_outfile)
        write_germline_vcf(vcf_list[target_id], germline_outfile)














    
    
            
        
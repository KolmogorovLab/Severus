#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from collections import defaultdict
from datetime import datetime
import sys
import os
import math




class vcf_format(object):
    __slots__ = ('chrom', 'pos', 'haplotype', 'ID', 'sv_type','alt', 'sv_len', 'qual', 'Filter', 'chr2', 'pos2', 'DR', 'DV', 'mut_type', 'hVaf', 'ins_seq', 'gen_type', 'cluster_id', 'ins_len', 'detailed_type', 'prec')
    def __init__(self, chrom, pos, haplotype, ID, sv_type, alt, sv_len, qual, Filter, chr2, pos2, DR, DV, mut_type, hVaf, gen_type, cluster_id, ins_len, detailed_type, prec):
        self.chrom = chrom
        self.pos = pos
        self.haplotype = haplotype
        self.ID = ID
        self.sv_type = sv_type
        self.alt = alt
        self.sv_len = sv_len
        self.qual = qual
        self.Filter = Filter
        self.chr2 = chr2
        self.pos2 = pos2
        self.DR = DR
        self.DV = DV
        self.hVaf = hVaf
        self.mut_type = mut_type
        self.ins_seq = None
        self.gen_type = gen_type
        self.cluster_id = cluster_id
        self.ins_len = ins_len
        self.detailed_type = detailed_type
        self.prec = prec
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
        GT = genotypes[g_list.index(min(g_list))]
        if GT == "0/1" and self.gen_type:
            GT = '1|0' if self.gen_type == 'HP1' else '0|1'
        g_list.sort()
        GQ = int(g_list[2])
        return GT, GQ
    def precision(self):
        return 'PRECISE' if self.prec else 'IMPRECISE' 
    def info(self):
        return f"{self.precision()};SVTYPE={self.sv_type};SVLEN={self.sv_len};CHR2={self.chr2};END={self.pos2};DETAILED_TYPE={self.detailed_type};INSLEN={self.ins_len};MAPQ={self.qual};SUPPREAD={self.DV};HVAF={self.hVAF()};CLUSTERID=severus_{self.cluster_id}"
    def sample(self):
        GT, GQ= self.call_genotype()
        return f"{GT}:{GQ}:{self.vaf():.2f}:{self.DR}:{self.DV}"
    def to_vcf(self):
        if self.ins_seq:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t<{self.ins_seq}>\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:GQ:VAF:DR:DV\t{self.sample()}\n"
        elif not self.sv_type == self.alt:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t{self.alt}\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:GQ:VAF:DR:DV\t{self.sample()}\n"
        else:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t<{self.sv_type}>\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:GQ:VAF:DR:DV\t{self.sample()}\n"
            
            
 
def get_sv_type(db):
    if db.bp_2.loose_end_id:
        db.sv_type = 'loose_end'
        return 'BND'
    if db.is_dup:
        return 'DUP'
    if db.bp_2.is_insertion or db.bp_1.is_insertion:
        return 'INS'
    if not db.bp_1.ref_id == db.bp_2.ref_id:
        return 'BND'
    if db.direction_1 == db.direction_2:
        return 'INV'
    if db.bp_1.dir_1 < 0:
        return 'INS'
    if db.bp_1.dir_1 > 0:
        return 'DEL'

def db_2_vcf(double_breaks, id_to_cc, no_ins):
    t = 0
    vcf_list = defaultdict(list)
    clusters = defaultdict(list) 
    for br in double_breaks:
        if br.bp_1.is_insertion:
            continue
        clusters[br.to_string()].append(br)
        
    for key, db_clust in clusters.items():
        cluster_id = id_to_cc[key]
        db_list = defaultdict(list)
        pass_list = [db.is_pass for db in db_clust]
        new_pass = True if 'PASS' in pass_list else False
        
        for db in db_clust:
            if new_pass:
                db.is_pass = 'PASS'
            db_list[db.genome_id].append(db)
            
        sv_type = get_sv_type(db_clust[0])
        if not sv_type:
            continue
        
        for db1 in db_list.values():
            if db1[0].genotype == 'hom':
                gen_type = ''
            else:
                hap_type = [db.haplotype_1 for db in db1]
                gen_type = 'HP1' if 1 in hap_type else 'HP2'
            ID = 'SEVERUS_' + sv_type + str(t)
            t += 1
            qual_list = []
            pass_list = [db.is_pass for db in db1]
            if 'PASS' in pass_list:
                sv_pass = 'PASS' 
            elif 'FAIL_LOWCOV_OTHER' in pass_list:
                sv_pass = 'FAIL_LOWCOV_OTHER' 
            else:
                sv_pass = db1[0].is_pass
            hVaf = [0.0, 0.0, 0.0]
            for db in db1:
                if db.is_pass == 'PASS':
                    qual_list += [db.bp_1.qual, db.bp_2.qual]
                hVaf[db.haplotype_1]  = db.hvaf
            if not qual_list:
                qual_list = [0]
            if sv_type == "BND" and not db.bp_2.loose_end_id:
                if db.bp_2.dir_1 == 1:
                    alt = '[' + db.bp_2.ref_id +':'+ str(db.bp_2.position) + '[N'
                else:
                    alt = 'N]' + db.bp_2.ref_id +':'+ str(db.bp_2.position) + '['
                vcf_list[db.genome_id].append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, alt, db.length, int(np.median(qual_list)), 
                                                         sv_pass, db.bp_2.ref_id, db.bp_2.position, db.DR, db.DV, db.mut_type, hVaf, gen_type, cluster_id,db.has_ins,db.sv_type, db.prec))#
                if db.bp_1.dir_1 == 1:
                    alt = '[' + db.bp_1.ref_id +':'+ str(db.bp_1.position) + '[N'
                else:
                    alt = 'N]' + db.bp_1.ref_id +':'+ str(db.bp_1.position) + '['
                    
                vcf_list[db.genome_id].append(vcf_format(db.bp_2.ref_id, db.bp_2.position, db.haplotype_1, ID, sv_type, alt, db.length, int(np.median(qual_list)), 
                                                         sv_pass, db.bp_1.ref_id, db.bp_1.position, db.DR, db.DV, db.mut_type, hVaf, gen_type, cluster_id,db.has_ins,db.sv_type, db.prec))#
            else:
            
                vcf_list[db.genome_id].append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, sv_type, db.length, int(np.median(qual_list)), 
                                                     sv_pass, db.bp_2.ref_id, db.bp_2.position, db.DR, db.DV, db.mut_type, hVaf, gen_type, cluster_id,db.has_ins,db.sv_type, db.prec))#
            
            if sv_type == 'INS':
                if not no_ins:
                    vcf_list[db.genome_id][-1].ins_seq = db.ins_seq
                else:
                    vcf_list[db.genome_id][-1].ins_seq = 'INS'
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

    outfile.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    outfile.write('##FILTER=<ID=FAIL_LOWSUPP,Description="Less number of support, but ok in other samples">\n')
    outfile.write('##FILTER=<ID=FAIL_MAP_CONS,Description="Majority of variant reads have unreliable mappability">\n')
    outfile.write('##FILTER=<ID=FAIL_CONN_CONS,Description="Majority of variant reads have unreliable connections">\n')
    outfile.write('##FILTER=<ID=FAIL_LOWCOV_OTHER,Description="Low variant coverage in other samples">\n')

    outfile.write("##INFO=<ID=HAPLOTYPE,Number=1,Type=String,Description=\"Haplotype of the SV\">\n")
    outfile.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    outfile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    outfile.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">\n")
    outfile.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV\">\n")
    outfile.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of supporting reads\">\n")
    outfile.write("##INFO=<ID=SUPPREAD,Number=1,Type=Integer,Description=\"Number of supporting reads\">\n")#
    outfile.write("##INFO=<ID=HVAF,Number=1,Type=String,Description=\"Haplotype specific variant Allele frequency (H0|H1|H2)\">\n")
    outfile.write("##INFO=<ID=CLUSTERID,Number=1,Type=String,Description=\"Cluster ID in breakpoint_graph\">\n")

    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    outfile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotyping quality\">\n")
    outfile.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of reference reads\">\n")
    outfile.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of variant reads\">\n")
    outfile.write("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    

def _sorted_breakpoints(bp_list):
    def chr_sort(x):
        stripped_chr = x[3:] if x.startswith("chr") else x
        return int(stripped_chr) if stripped_chr.isnumeric() else sum(map(ord, stripped_chr))

    return sorted(bp_list, key=lambda x: (chr_sort(x.chrom), x.pos))


def write_somatic_vcf(vcf_list, outfile):
    for db in _sorted_breakpoints(vcf_list):
        if db.mut_type == 'somatic':
            outfile.write(db.to_vcf())
    outfile.close()

    
def write_germline_vcf(vcf_list, outfile):
    for db in _sorted_breakpoints(vcf_list):
        outfile.write(db.to_vcf())
    outfile.close()
    
    
def write_to_vcf(double_breaks, all_ids, id_to_cc, outpath, out_key, ref_lengths, only_somatic, no_ins):
    vcf_list = db_2_vcf(double_breaks, id_to_cc, no_ins)
    key = 'somatic_' if out_key == 'somatic' else 'all_'
    for target_id in all_ids:
        germline_outfile = open(os.path.join(outpath, 'severus_' + key + target_id.replace('.bam' , '') + ".vcf"), "w")
        write_vcf_header(ref_lengths, germline_outfile)
        write_germline_vcf(vcf_list[target_id], germline_outfile)
    
            
    


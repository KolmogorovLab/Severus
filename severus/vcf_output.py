#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from collections import defaultdict
from datetime import datetime
import sys
import os
import math




class vcf_format(object):
    __slots__ = ('chrom', 'pos', 'haplotype', 'ID', 'sv_type','alt', 'sv_len', 'qual', 'Filter', 'chr2', 'pos2','mut_type',
                 'ins_seq', 'ins_len_seq', 'cluster_id','ins_len', 'detailed_type', 'prec', 'phaseset_id', 'strands', 'sample', 'gen_type')
    def __init__(self, chrom, pos, haplotype, ID, sv_type, alt, sv_len, qual, Filter, chr2, pos2, mut_type, cluster_id, 
                 ins_len, ins_len_seq, detailed_type, prec, phaseset_id, strands, sample, gen_type):
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
        self.mut_type = mut_type
        self.ins_seq = None
        self.cluster_id = cluster_id
        self.ins_len = ins_len
        self.ins_len_seq = ''
        self.detailed_type = detailed_type
        self.prec = prec
        self.phaseset_id = phaseset_id
        self.strands = strands
        self.sample = sample
        self.gen_type = gen_type
     
    def precision(self):
        return 'PRECISE' if self.prec else 'IMPRECISE'
    
    def strands_vcf(self):
        if self.sv_type == 'INS':
            return ''
        else:
            return f"STRANDS={self.strands[0]}{self.strands[1]};"
        
    def has_ins(self):
        if self.ins_len:
            return f"INSLEN={self.ins_len};INSSEQ={self.ins_len_seq};"
        else:
            return ""
        
    def detailedtype(self):
        if self.detailed_type:
            return f"DETAILED_TYPE={self.detailed_type};"
        else:
            return ""
        
    def info(self):
        if self.gen_type and self.phaseset_id:
            phase_id = str(self.phaseset_id[0]) if self.phaseset_id[0] == self.phaseset_id[1] else '{0}|{1}'.format(self.phaseset_id[0], self.phaseset_id[1])
            return f"{self.precision()};SVTYPE={self.sv_type};SVLEN={self.sv_len};CHR2={self.chr2};END={self.pos2};{self.strands_vcf()}{self.detailedtype()}{self.has_ins()}MAPQ={self.qual};PHASESET_ID={phase_id};CLUSTERID=severus_{self.cluster_id}"
        else:
            return f"{self.precision()};SVTYPE={self.sv_type};SVLEN={self.sv_len};CHR2={self.chr2};END={self.pos2};{self.strands_vcf()}{self.detailedtype()}{self.has_ins()}MAPQ={self.qual};CLUSTERID=severus_{self.cluster_id}"
    
    def to_vcf(self):
        if self.ins_seq:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t{self.ins_seq}\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:GQ:VAF:hVAF:DR:DV\t{self.sample}\n"
        elif not self.sv_type == self.alt:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t{self.alt}\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:GQ:VAF:hVAF:DR:DV\t{self.sample}\n"
        else:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t<{self.sv_type}>\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:GQ:VAF:hVAF:DR:DV\t{self.sample}\n"
      
        
class vcf_sample(object):
    __slots__ = ('DR', 'DV', 'vaf', 'hVaf','gen_type', 'GT', 'GQ')
    def __init__(self, DR, DV, vaf, hVaf, gen_type):
        self.DR = DR
        self.DV = DV
        self.hVaf = hVaf
        self.vaf = vaf
        self.gen_type = gen_type
        self.GT = None
        self.GQ = None
        
    def hVAF(self):
        return '{0:.2f},{1:.2f},{2:.2f}'.format(self.hVaf[0], self.hVaf[1], self.hVaf[2])
    
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
        self.GT = GT
        self.GQ = GQ
    
    def sample(self):
        if self.DV == self.DR == 0:
            return "0:0:0:0,0,0:0:0"
        self.call_genotype()
        return f"{self.GT}:{self.GQ}:{self.vaf:.2f}:{self.hVAF()}:{self.DR}:{self.DV}"
    
    
def db_2_vcf(double_breaks, no_ins, sample_ids):
    vcf_list = []
    clusters = defaultdict(list) 
    for br in double_breaks:
        if br.bp_1.is_insertion:
            continue
        clusters[br.to_string()].append(br)
        
    for key, db_clust in clusters.items():
        db_list = defaultdict(list)
        
        pass_list = [db.is_pass for db in db_clust]
        new_pass = True if 'PASS' in pass_list else False
        
        for db in db_clust:
            if new_pass:
                db.is_pass = 'PASS'
            db_list[db.genome_id].append(db)
            
        sv_type = db.vcf_sv_type
        
        pass_list = [db.is_pass for db in db_clust]
        if 'PASS' in pass_list:
            sv_pass = 'PASS' 
        elif 'FAIL_LOWCOV_OTHER' in pass_list:
            sv_pass = 'FAIL_LOWCOV_OTHER'
        else:
            sv_pass = db[0].is_pass
        
        dir1 = '-' if db.direction_1 == -1 else '+'
        dir2 = '-' if db.direction_2 == -1 else '+'
        strands = (dir1, dir2)
        
        ID = db.vcf_id
        
        if db.genotype == 'hom':
            gen_type = ''
            phaseset = ''
        else:
            hap_type = [db.haplotype_1 for db in db_clust]
            gen_type = 'HP1' if 1 in hap_type else 'HP2'
            phaseset = db.phaseset_id
        
        sample_list = defaultdict(list)
        for sample_id in sample_ids:
            sample_list[sample_id] = vcf_sample(0,0,0,0,0)
        for db1 in db_list.values():
            hVaf = [0.0, 0.0, 0.0]
            for db in db1:
                hVaf[db.haplotype_1]  = db.hvaf
            sample_list[db.genome_id] = vcf_sample(db.DR, db.DV, db.vaf, hVaf, gen_type)
          
        sample = '\t'.join([s.sample() for s in sample_list.values()])
        phased = None
        GT = [s.GT for s in sample_list.values()]
        if '1|0' in GT or '0|1' in GT:
            phased= True
            
        if sv_type == "BND" and not db.bp_2.loose_end_id:
            if db.bp_2.dir_1 == 1 and db.bp_1.dir_1 == 1:
                alt1 = 'N]' + db.bp_2.ref_id +':'+ str(db.bp_2.position) + ']'
                alt2 = 'N]' + db.bp_1.ref_id +':'+ str(db.bp_1.position) + ']'
            elif db.bp_2.dir_1 == -1 and db.bp_1.dir_1 == -1:
                alt1 = '[' + db.bp_2.ref_id +':'+ str(db.bp_2.position) + '[N'
                alt2 = '[' + db.bp_1.ref_id +':'+ str(db.bp_1.position) + '[N'
            elif db.bp_2.dir_1 == -1 and db.bp_1.dir_1 == 1:
                alt1 = ']' + db.bp_2.ref_id +':'+ str(db.bp_2.position) + ']N'
                alt2 = 'N[' + db.bp_1.ref_id +':'+ str(db.bp_1.position) + '['
            else:
                alt1 = 'N[' + db.bp_2.ref_id +':'+ str(db.bp_2.position) + '['
                alt2 = ']' + db.bp_1.ref_id +':'+ str(db.bp_1.position) + ']N'
            
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, alt1, db.length, db.vcf_qual, 
                                                     sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type,db.cluster_id,
                                                     db.has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, phased))#
            vcf_list.append(vcf_format(db.bp_2.ref_id, db.bp_2.position, db.haplotype_1, ID, sv_type, alt2, db.length, db.vcf_qual, 
                                                     sv_pass, db.bp_1.ref_id, db.bp_1.position, db.mut_type, db.cluster_id,
                                                     db.has_ins, db.ins_seq,db.sv_type, db.prec, phaseset, strands,sample, phased))
        else:
        
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, sv_type, db.length, db.vcf_qual, 
                                                 sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type, db.cluster_id,
                                                 db.has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, phased))#
        
        if sv_type == 'INS':
            if not no_ins:
                vcf_list[-1].ins_seq = db.ins_seq
            else:
                vcf_list[-1].ins_seq = 'INS'
                
    return vcf_list


          
def write_vcf_header(ref_lengths, outfile, sample_list):
    sample = '\t'.join(sample_list)
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write('##source=Severusv0.1.2\n')
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

    outfile.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"SV with precise breakpoints coordinates and length\">\n")
    outfile.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"SV with imprecise breakpoints coordinates and length\">\n")
    outfile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    outfile.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    outfile.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">\n")
    outfile.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV\">\n")
    outfile.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strandedness\">\n")
    outfile.write("##INFO=<ID=DETAILED_TYPE,Number=1,Type=Integer,Description=\"Detailed type of the SV\">\n")
    outfile.write("##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Length of the unmapped sequence between breakpoint\">\n")
    outfile.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of supporting reads\">\n")
    outfile.write("##INFO=<ID=PHASESET_ID,Number=1,Type=Integer,Description=\"Matching phaseset ID for phased SVs\">\n")
    outfile.write("##INFO=<ID=CLUSTERID,Number=1,Type=String,Description=\"Cluster ID in breakpoint_graph\">\n")
    outfile.write("##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Insertion sequence between breakpoints\">\n")

    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    outfile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotyping quality\">\n")
    outfile.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of reference reads\">\n")
    outfile.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of variant reads\">\n")
    outfile.write("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n")
    outfile.write("##FORMAT=<ID=hVAF,Number=3,Type=Float,Description=\"Haplotype specific variant Allele frequency (H0,H1,H2)\">\n")
    outfile.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
    

def _sorted_breakpoints(bp_list):
    def chr_sort(x):
        stripped_chr = x[3:] if x.startswith("chr") else x
        return int(stripped_chr) if stripped_chr.isnumeric() else sum(map(ord, stripped_chr))

    return sorted(bp_list, key=lambda x: (chr_sort(x.chrom), x.pos))

    
def write_germline_vcf(vcf_list, outfile):
    for db in _sorted_breakpoints(vcf_list):
        outfile.write(db.to_vcf())
    outfile.close()
    
    
def write_to_vcf(double_breaks, all_ids, outpath, out_key, ref_lengths, only_somatic, no_ins):
    vcf_list = db_2_vcf(double_breaks, no_ins, all_ids)
    key = 'somatic' if out_key == 'somatic' else 'all'
    sample_ids = [target_id.replace('.bam' , '') for target_id in all_ids]
    
    germline_outfile = open(os.path.join(outpath, 'severus_' + key + ".vcf"), "w")
    write_vcf_header(ref_lengths, germline_outfile, sample_ids)
    write_germline_vcf(vcf_list, germline_outfile)
    
            
    


#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from collections import defaultdict
from datetime import datetime
import sys
import os
import math




class vcf_format(object):
    __slots__ = ('chrom', 'pos', 'haplotype', 'ID', 'sv_type','alt', 'sv_len', 'qual', 'Filter', 'chr2', 'pos2','mut_type', 'tra_pos',
                 'ins_seq', 'ins_len_seq', 'cluster_id','ins_len', 'detailed_type', 'prec', 'phaseset_id', 'strands', 'sample','HP', 'mate_id', 'vntr')
    def __init__(self, chrom, pos, haplotype, ID, sv_type, alt, sv_len, qual, Filter, chr2, pos2, mut_type, cluster_id, 
                 ins_len, ins_len_seq, detailed_type, prec, phaseset_id, strands, sample, HP, mate_id, vntr,tra_pos):
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
        self.ins_len_seq = ins_len_seq
        self.detailed_type = detailed_type
        self.prec = prec
        self.phaseset_id = phaseset_id
        self.strands = strands
        self.sample = sample
        self.HP = HP
        self.mate_id = mate_id
        self.vntr = vntr
        self.tra_pos = tra_pos
     
    def precision(self):
        return 'PRECISE' if self.prec else 'IMPRECISE'
    
    def strands_vcf(self):
        if self.sv_type == 'INS':
            return ''
        if len(self.strands) == 1:
            return f"STRANDS={self.strands[0]};"
        else:
            return f"STRANDS={self.strands[0]}{self.strands[1]};"
        
    def has_ins(self):
        if self.ins_len:
            return f"INSLEN={self.ins_len};INSSEQ={self.ins_len_seq};"
        else:
            return ""
    
    def add_vntr(self):
        if self.vntr:
            return "INSIDE_VNTR=TRUE;"
        else:
            return ''
    def detailedtype(self):
        if self.detailed_type:
            return f"DETAILED_TYPE={self.detailed_type};"
        else:
            return ""
    def trapos(self):
        if self.tra_pos:
            return f"ALINGED_POS={self.tra_pos};"
        else:
            return ""
    def clusterid(self):
        if self.cluster_id and 'severus' in self.cluster_id:
            return f";CLUSTERID={self.cluster_id}"
        else:
            return ""
    def end_pos(self):
        if self.sv_type == 'BND':
            if self.mate_id:
                return f"CHR2={self.chr2};END={self.pos2};MATE_ID={self.mate_id};"
            else:
                return ""
        if self.ins_seq:
            return ""
        else:
            return f"CHR2={self.chr2};END={self.pos2};"
        
    def svlen(self):
        if self.sv_len:
            return f"SVLEN={self.sv_len};"
        else:
            return ""
            
    def HP_inf(self):
        if self.HP and self.phaseset_id:
            phase_id = str(self.phaseset_id[0]) if self.phaseset_id[0] == self.phaseset_id[1] else '{0}|{1}'.format(self.phaseset_id[0], self.phaseset_id[1])
            return f";PHASESETID={phase_id};HP={self.HP}"
        else:
            return ""
    def info(self):
        return f"{self.precision()};SVTYPE={self.sv_type};{self.svlen()}{self.end_pos()}{self.trapos()}{self.strands_vcf()}{self.add_vntr()}{self.detailedtype()}{self.has_ins()}MAPQ={self.qual}{self.HP_inf()}{self.clusterid()}"
    
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
        
        #if self.gen_type:
        #    GT = '1|0' if self.gen_type == 'HP1' else '0|1'
            
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
        if db.is_single:
            strands = (dir1)
        
        
        ID = db.vcf_id
        vntr = db.vntr
        
        if db.genotype == 'hom':
            gen_type1 = ''
            phaseset = ''
            gen_type2 = ''
        else:
            hap_type = [db.haplotype_1 for db in db_clust]
            if sum(hap_type) == 0 or sum(hap_type) == 3:
                gen_type1 = ''
            else:
                gen_type1 = '1' if 1 in hap_type else '2'
            
            hap_type = [db.haplotype_2 for db in db_clust]
            if sum(hap_type) == 0 or sum(hap_type) == 3:
                gen_type2 = ''
            else: 
                gen_type2 = '1' if 1 in hap_type else '2'
            phaseset = db.phaseset_id
        
        sample_list = defaultdict(list)
        sample_list2 = defaultdict(list)
        for sample_id in sample_ids:
            sample_list[sample_id] = vcf_sample(0,0,0,0,0)
            sample_list2[sample_id] = vcf_sample(0,0,0,0,0)
        for db1 in db_list.values():
            hVaf1 = [0.0, 0.0, 0.0]
            hVaf2 = [0.0, 0.0, 0.0]
            for db in db1:
                hVaf1[db.haplotype_1]  = db.hvaf
                hVaf2[db.haplotype_1]  = db.hvaf
            sample_list[db.genome_id] = vcf_sample(db.DR, db.DV, db.vaf, hVaf1, gen_type1)
            sample_list2[db.genome_id] = vcf_sample(db.DR, db.DV, db.vaf, hVaf2, gen_type2)
            
          
        sample = '\t'.join([s.sample() for s in sample_list.values()])
        sample2 = '\t'.join([s.sample() for s in sample_list2.values()])
            
        if sv_type == "BND" and not db.is_single:
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
                
            ID1 = ID + '_1'
            ID2 = ID + '_2'
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID1, sv_type, alt1, db.length, db.vcf_qual, 
                                                     sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type,db.cluster_id,
                                                     db.has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, gen_type1, ID2,vntr, db.tra_pos))#
            vcf_list.append(vcf_format(db.bp_2.ref_id, db.bp_2.position, db.haplotype_1, ID2, sv_type, alt2, db.length, db.vcf_qual, 
                                                     sv_pass, db.bp_1.ref_id, db.bp_1.position, db.mut_type, db.cluster_id,
                                                     db.has_ins, db.ins_seq,db.sv_type, db.prec, phaseset, strands,sample2, gen_type2, ID1,vntr, db.tra_pos))
        elif db.is_single:
            alt= 'N'
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, alt, db.length, db.vcf_qual, 
                                                 sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type, db.cluster_id,
                                                 db.has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, gen_type1, None,vntr, db.tra_pos))
            
        else:
        
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, db.haplotype_1, ID, sv_type, sv_type, db.length, db.vcf_qual, 
                                                 sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type, db.cluster_id,
                                                 db.has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, gen_type1, None,vntr, db.tra_pos))#
        
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
    outfile.write("##INFO=<ID=DETAILED_TYPE,Number=1,Type=String,Description=\"Detailed type of the SV\">\n")
    outfile.write("##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Length of the unmapped sequence between breakpoint\">\n")
    outfile.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of supporting reads\">\n")
    outfile.write("##INFO=<ID=PHASESETID,Number=1,Type=String,Description=\"Matching phaseset ID for phased SVs\">\n")
    outfile.write("##INFO=<ID=HP,Number=1,Type=Integer,Description=\"Matching haplotype ID for phased SVs\">\n")
    outfile.write("##INFO=<ID=CLUSTERID,Number=1,Type=String,Description=\"Cluster ID in breakpoint_graph\">\n")
    outfile.write("##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Insertion sequence between breakpoints\">\n")
    outfile.write("##INFO=<ID=MATE_ID,Number=1,Type=String,Description=\"MATE ID for breakends\">\n")
    outfile.write("##INFO=<ID=INSIDE_VNTR,Number=1,Type=String,Description=\"True if an indel is inside a VNTR\">\n")
    outfile.write("##INFO=<ID=ALINGED_POS,Number=1,Type=String,Description=\"Position in the reference\">\n")

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
    
    
def write_to_vcf(double_breaks, all_ids, outpath, out_key, ref_lengths, no_ins):
    vcf_list = db_2_vcf(double_breaks, no_ins, all_ids)
    key = 'somatic' if out_key == 'somatic' else 'all'
    sample_ids = [target_id.replace('.bam' , '') for target_id in all_ids]
    
    germline_outfile = open(os.path.join(outpath, 'severus_' + key + ".vcf"), "w")
    write_vcf_header(ref_lengths, germline_outfile, sample_ids)
    write_germline_vcf(vcf_list, germline_outfile)
    
            
    


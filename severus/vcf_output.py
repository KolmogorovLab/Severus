#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from severus.__version__ import __version__
from collections import defaultdict
from datetime import datetime
import numpy as np
import sys
import os
import math

class vcf_format(object):
    __slots__ = ('chrom', 'pos', 'haplotypes','haplotype', 'ID', 'sv_type','alt', 'sv_len', 'qual', 
                 'Filter', 'chr2', 'pos2','mut_type', 'tra_pos','ins_seq', 'ins_len_seq', 
                 'cluster_id','ins_len', 'detailed_type', 'prec', 'phaseset_id', 'strands', 
                 'sample','HP', 'mate_id', 'vntr', 'low_cov', 'whitelist', 'dvls', 'drls', 'bnd_type',
                 'wl_rescue')
    def __init__(self, chrom, pos, haplotypes, haplotype, ID, sv_type, alt, sv_len, qual, 
                 Filter, chr2, pos2, mut_type, cluster_id, ins_len, ins_len_seq, detailed_type, 
                 prec, phaseset_id, strands, sample, HP, mate_id, vntr,tra_pos, low_cov, whitelist,
                 dvls, drls,bnd_type, wl_rescue=None):
        self.chrom = chrom
        self.pos = pos
        self.haplotypes = haplotypes
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
        self.low_cov = low_cov
        self.whitelist = whitelist
        self.dvls = dvls
        self.drls = drls
        self.bnd_type = bnd_type
        self.wl_rescue = wl_rescue
     
    def precision(self):
        return 'PRECISE' if self.prec else 'IMPRECISE'
    def inwhitelist(self):
        if self.whitelist:
            return ';INSIDE_WHITELIST=TRUE'
        else:
            return ''
    def wlrescue(self):
        if self.wl_rescue:
            return ';WL_RESCUE=' + self.wl_rescue
        else:
            return ''
    
    def strands_vcf(self):
        if self.sv_type == 'INS':
            return ''
        if len(self.strands) == 1:
            return f"STRANDS={self.strands[0]};"
        else:
            return f"STRANDS={self.strands[0]}{self.strands[1]};"
        
    def has_ins(self):
        if self.sv_type == 'INS':
            return ""
        if self.ins_len and self.ins_len_seq:
            return f"INSLEN={self.ins_len};INSSEQ={self.ins_len_seq};"
        elif self.ins_len:
            return f"INSLEN={self.ins_len};"
        else:
            return ""
    
    def add_vntr(self):
        if self.vntr:
            return "INSIDE_VNTR=TRUE;"
        else:
            return ''
    def detailedtype(self):
        if self.detailed_type and not self.detailed_type == 'DEL':
            return f"DETAILED_TYPE={self.detailed_type};"
        else:
            return ""
    def trapos(self):
        if self.tra_pos:
            return f"ALIGNED_POS={self.tra_pos};"
        else:
            return ""
    def clusterid(self):
        if self.cluster_id and 'severus' in self.cluster_id:
            return f";CLUSTERID={self.cluster_id}"
        else:
            return ""
    def end_pos(self):
        if self.alt == '.N':
            return ''
        if self.sv_type == 'BND':
            if self.mate_id:
                return f"MATE_ID={self.mate_id};"
            else:
                return ""
        elif self.sv_type == 'DEL' or self.sv_type == 'DUP' or self.sv_type == 'INV':
            return  f"END={self.pos2};"
        else:
            return ""
        
    def svlen(self):
        if self.sv_len > 0:
            return f"SVLEN={self.sv_len};"
        else:
            return ""
        
    def HP_inf(self):
        if self.haplotypes and self.phaseset_id:
            phase_id = '{0}|{1}'.format(self.phaseset_id[0], self.phaseset_id[1])
            hps = '{0}|{1}'.format(self.haplotypes[0], self.haplotypes[1])
            return f";PHASESETID={phase_id};HP={hps}"
        else:
            return ""
    def to_dvls(self):
        return ";SUPP_READS=" + ':'.join(self.dvls[0]) + ":"+ ':'.join(self.dvls[1]) +";REF_READS=" + ':'.join(self.drls[0]) + ":"+ ':'.join(self.drls[1])

    def low_covls(self):
        return 'LOW_COV_IN=' + ','.join(self.low_cov) + ';' if self.low_cov else ''
    
    def bndtype(self):
        return 'BND_TYPE=' + self.bnd_type + ';' if self.bnd_type else ''
    
    def info(self):
        return f"{self.precision()};SVTYPE={self.sv_type};{self.svlen()}{self.end_pos()}{self.trapos()}{self.bndtype()}{self.strands_vcf()}{self.low_covls()}{self.add_vntr()}{self.detailedtype()}{self.has_ins()}MAPQ={self.qual}{self.HP_inf()}{self.to_dvls()}{self.clusterid()}{self.inwhitelist()}{self.wlrescue()}"
    
    def to_vcf(self):
        if self.sv_type == 'INS':
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t{self.ins_seq}\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:VAF:hVAF:DR:DV\t{self.sample}\n"
        elif not self.sv_type == self.alt:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t{self.alt}\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:VAF:hVAF:DR:DV\t{self.sample}\n"
        else:
            return f"{self.chrom}\t{self.pos}\t{self.ID}\tN\t<{self.sv_type}>\t{self.qual}\t{self.Filter}\t{self.info()}\tGT:VAF:hVAF:DR:DV\t{self.sample}\n"
      
        
class vcf_sample(object):
    __slots__ = ('DR', 'DV', 'vaf', 'hVaf','gen_type', 'GT', 'GQ', 'genotype')
    def __init__(self, DR, DV, vaf, hVaf, gen_type, genotype):
        self.DR = DR
        self.DV = DV
        self.hVaf = hVaf
        self.vaf = vaf
        self.gen_type = gen_type
        self.GT = None
        self.GQ = None
        self.genotype = genotype
        
    def hVAF(self):
        return '{0:.2f},{1:.2f},{2:.2f}'.format(self.hVaf[0], self.hVaf[1], self.hVaf[2])
        
    
    def call_genotype(
        self,
        err: float = 0.01,
        priors=(0.90, 0.09, 0.01),
        max_read: int = 100,
        pl_cap: int = 999,
        force_phase_if_het: bool = True,
    ):
        # --- 1) Get counts, cap to avoid extreme numbers dominating & reduce numeric issues
        DV = int(getattr(self, "DV", 0) or 0)
        DR = int(getattr(self, "DR", 0) or 0)

        if DV < 0 or DR < 0:
            # Defensive: invalid counts
            self.GT, self.GQ, self.PL = "./.", 0, [pl_cap, pl_cap, pl_cap]
            return self.GT, self.GQ, self.PL

        tot = DV + DR
        if tot == 0:
            # No evidence: typically no-call or 0/0 low confidence; choose your preference
            self.GT, self.GQ, self.PL = "./.", 0, [0, 0, 0]
            return self.GT, self.GQ, self.PL

        if tot > max_read:
            DV = int(round(max_read * DV / tot))
            DR = max_read - DV
            tot = max_read

        # --- 2) Log-likelihoods: P(data | GT) * prior(GT)
        p00, p01, p11 = priors
        # guard against priors=0
        p00 = max(p00, 1e-300)
        p01 = max(p01, 1e-300)
        p11 = max(p11, 1e-300)

        # guard against err extremes
        err = min(max(err, 1e-12), 1 - 1e-12)

        # 0/0: variant reads are errors
        log_hom0 = DV * math.log(err) + DR * math.log(1 - err) + math.log(p00)
        # 0/1: symmetric 50/50 model
        log_het  = tot * math.log(0.5) + math.log(p01)
        # 1/1: ref reads are errors
        log_hom1 = DV * math.log(1 - err) + DR * math.log(err) + math.log(p11)

        logs = [log_hom0, log_het, log_hom1]
        m = max(logs)
        rel = [math.exp(x - m) for x in logs]  # in (0,1]
        pls = []
        for p in rel:
            if p <= 0.0:
                pls.append(pl_cap)
            else:
                pls.append(-10.0 * math.log10(p))

        min_pl = min(pls)
        pls = [min(pl_cap, pl - min_pl) for pl in pls]

        PL = [int(round(x)) for x in pls]

        gt_map = ["0/0", "0/1", "1/1"]
        best_idx = min(range(3), key=lambda i: PL[i])
        GT = gt_map[best_idx]

        sorted_pl = sorted(PL)
        GQ = int(sorted_pl[1])

        if force_phase_if_het:
            gt_phase = getattr(self, "gen_type", None)
            if gt_phase == "HP1":
                GT = "1|0"
            elif gt_phase == "HP2":
                GT = "0|1"

        self.GT, self.GQ = GT, GQ
    
    def sample(self):
        if self.DV == 0:
            return f"./.:0:0,0,0:{self.DR}:0"
        if self.genotype:
            self.call_genotype()
        else:
            self.GT = '0/1'
        return f"{self.GT}:{self.vaf:.2f}:{self.hVAF()}:{self.DR}:{self.DV}"
    
    
def db_2_vcf(double_breaks, no_ins, sample_ids, multisample, junction_vcf, germline_genotype, ignore_hp):
    vcf_list = []
    clusters = defaultdict(list)
    bnd_types = {'-+': 'DUP_LIKE', '+-': 'DEL_LIKE', '++': 'INV_LIKE', '--': 'INV_LIKE'}
        
    for br in double_breaks:
        if br.bp_1.is_insertion:
            continue
        clusters[br.to_string()].append(br)
        
    for key, db_clust in clusters.items():
        db_list = defaultdict(list)
        db_clust.sort(key=lambda k:k.supp, reverse=True)
        pass_list = [db.is_pass for db in db_clust]
        new_pass = True if 'PASS' in pass_list else False
        bnd_type = ''
        
        for db in db_clust:
            if new_pass:
                db.is_pass = 'PASS'
            if db.ins_seq == '<DUP>':
                db.ins_seq = ''
                db.has_ins = ''
            db_list[db.genome_id].append(db)
            
        sv_type = db.vcf_sv_type
        
        pass_list = [db.is_pass for db in db_clust]
        if 'PASS' in pass_list:
            sv_pass = 'PASS' 
        elif 'FAIL_LOWCOV_OTHER' in pass_list:
            sv_pass = 'FAIL_LOWCOV_OTHER'
        else:
            sv_pass = db_clust[0].is_pass
        if db.whitelist:
            qual_list = []
            for db in db_clust:
                qual_list += [db.bp_1.qual, db.bp_2.qual]
            vcf_qual = np.median(qual_list)
            for db in db_clust:
                db.vcf_qual = vcf_qual
        low_cov = None
        if multisample:
            low_cov = set()
            MIN_SPAN = 5
            db = db_clust[0]
            hpls1 = list(set([db.haplotype_1 for db in db_clust]))
            hpls2 = list(set([db.haplotype_1 for db in db_clust]))
            for genome_id, spann_ls in db.bp_1.spanning_reads.items():
                spann = 0
                for hp in hpls1:
                    spann += spann_ls[hp]
                if spann < MIN_SPAN:
                    low_cov.add(genome_id)
                    
            for genome_id, spann_ls in db.bp_2.spanning_reads.items():
                spann = 0
                for hp in hpls2:
                    spann += spann_ls[hp]
                if spann < MIN_SPAN:
                    low_cov.add(genome_id)
            
            low_cov = list(low_cov)
        
        dir1 = '-' if db.direction_1 == -1 else '+'
        dir2 = '-' if db.direction_2 == -1 else '+'
        strands = (dir1, dir2)
        strands2 = (dir2, dir1)
        if db.is_single:
            strands = (dir1)
        
        ID = db.vcf_id
    
        vntr = True if db.vntr else None
        
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
        haplotype = (db.haplotype_1, db.haplotype_2)
        if ignore_hp:
            phaseset = ''
            haplotype = None
        
        drls = [[str(s) for s in db.bp_1.spanning_reads[db.genome_id][:3]], [str(s) for s in db.bp_2.spanning_reads[db.genome_id][:3]]]
        dvls = [[str(s) for s in db.dvls[0]], [str(s) for s in db.dvls[1]]]
        
        sample_list = defaultdict(list)
        sample_list2 = defaultdict(list)
        for sample_id in sample_ids:
            sample_list[sample_id] = vcf_sample(sum(db.bp_1.spanning_reads[sample_id]),0,0,0,0, germline_genotype)
            sample_list2[sample_id] = vcf_sample(sum(db.bp_2.spanning_reads[sample_id]),0,0,0,0, germline_genotype)
        for db1 in db_list.values():
            hVaf1 = [0.0, 0.0, 0.0]
            hVaf2 = [0.0, 0.0, 0.0]
            for db in db1:
                hVaf1[db.haplotype_1]  = db.hVaf
                hVaf2[db.haplotype_1]  = db.hVaf
            sample_list[db.genome_id] = vcf_sample(db.DR, db.DV, db.vaf, hVaf1, gen_type1, germline_genotype)
            sample_list2[db.genome_id] = vcf_sample(db.DR, db.DV, db.vaf, hVaf2, gen_type2, germline_genotype)
            
          
        sample = '\t'.join([s.sample() for s in sample_list.values()])
        sample2 = '\t'.join([s.sample() for s in sample_list2.values()])
        
        has_ins = 0 if not db.ins_seq else len(db.ins_seq)
        if (sv_type == "BND"  and not db.is_single) or (junction_vcf and not sv_type == 'INS'):
            b1 = f"{db.bp_1.ref_id}:{db.bp_1.position}"
            b2 = f"{db.bp_2.ref_id}:{db.bp_2.position}"
        
            s1 = db.bp_1.dir_1  # local sign at bp1 (+1 or -1)
            s2 = db.bp_2.dir_1  # local sign at bp2 (+1 or -1)
            
            bnd_type = 'INTER_CHR' if not db.bp_1.ref_id == db.bp_2.ref_id else bnd_types[dir1+dir2]
        
            if s1 == 1 and s2 == 1:          # +/+
                alt1 = f"N]{b2}]"
                alt2 = f"N]{b1}]"
            elif s1 == 1 and s2 == -1:       # +/-
                alt1 = f"N[{b2}["
                alt2 = f"]{b1}]N"
            elif s1 == -1 and s2 == 1:       # -/+
                alt1 = f"]{b2}]N"
                alt2 = f"N[{b1}["
            else:                             # -/-
                alt1 = f"[{b2}[N"
                alt2 = f"[{b1}[N"
                
            ID1 = ID + '_1'
            ID2 = ID + '_2'
            haplotype2 = (db.haplotype_2, db.haplotype_1)
            dvls2 = [dvls[1], dvls[0]]
            drls2 = [drls[1], drls[0]]
            if sv_type == 'INV':
                db.sv_type = 'reciprocal_inversion'
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, haplotype, db.haplotype_1, ID1, sv_type, alt1, db.length, db.vcf_qual, 
                                                     sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type,db.cluster_id,
                                                     has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, gen_type1, ID2,vntr, 
                                                     db.tra_pos, low_cov, db.whitelist, dvls, drls, bnd_type, db.wl_rescue))#
            vcf_list.append(vcf_format(db.bp_2.ref_id, db.bp_2.position, haplotype2, db.haplotype_2, ID2, sv_type, alt2, db.length, db.vcf_qual, 
                                                     sv_pass, db.bp_1.ref_id, db.bp_1.position, db.mut_type, db.cluster_id,
                                                     has_ins, db.ins_seq,db.sv_type, db.prec, phaseset, strands2,sample2, gen_type2, ID1,vntr, 
                                                     db.tra_pos, low_cov,db.whitelist, dvls2, drls2, bnd_type, db.wl_rescue))
        elif db.is_single:
            alt= '.N' if db.direction_1 == -1 else 'N.'
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, haplotype, db.haplotype_1, ID, 'sBND', alt, db.length, db.vcf_qual, 
                                                 sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type, db.cluster_id,
                                                 has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, gen_type1, None,vntr, db.bp_1.pos2, 
                                                 low_cov,db.whitelist, dvls, drls, bnd_type, db.wl_rescue))
            vcf_list[-1].pos2 = db.bp_1.pos2
        elif sv_type == 'INV':
            if db.direction_1 == -1:
                continue
            pos1 , pos2 = int(min(db.sv_type)), int(max(db.sv_type))
            vcf_list.append(vcf_format(db.bp_1.ref_id, pos1, haplotype, db.haplotype_1, ID, sv_type, sv_type, pos2-pos1, db.vcf_qual, 
                                                 sv_pass, db.bp_2.ref_id, pos2, db.mut_type, db.cluster_id,
                                                 has_ins, db.ins_seq, 'reciprocal_inversion', db.prec, phaseset, strands,sample, gen_type1, 
                                                 None,vntr, None, low_cov,db.whitelist, dvls, drls, bnd_type, db.wl_rescue))
        else:
            if not sv_type == 'INS' and db.bp_1.ref_id == db.bp_2.ref_id:
                db.length = abs(db.bp_2.position - db.bp_1.position)
        
            vcf_list.append(vcf_format(db.bp_1.ref_id, db.bp_1.position, haplotype, db.haplotype_1, ID, sv_type, sv_type, db.length, db.vcf_qual, 
                                                 sv_pass, db.bp_2.ref_id, db.bp_2.position, db.mut_type, db.cluster_id,
                                                 has_ins, db.ins_seq, db.sv_type, db.prec, phaseset, strands,sample, gen_type1, None,vntr, 
                                                 db.tra_pos, low_cov,db.whitelist, dvls, drls, bnd_type, db.wl_rescue))#
        
        if sv_type == 'INS':
            if not no_ins:
                vcf_list[-1].ins_seq = db.ins_seq
            else:
                vcf_list[-1].ins_seq = 'INS'
                
    return vcf_list


          
def write_vcf_header(ref_lengths, outfile, sample_list):

    sample = '\t'.join(sample_list)
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write('##source=Severus_v'+ __version__ + '\n')
    outfile.write('##CommandLine= '+ " ".join(sys.argv[1:]) +'\n')
    filedate = str(datetime.now()).split(' ')[0]
    outfile.write('##fileDate='+filedate+'\n')#
    for chr_id, chr_len in ref_lengths.items():
        outfile.write("##contig=<ID={0},length={1}>\n".format(chr_id, chr_len))#
    outfile.write('##ALT=<ID=DEL,Description="Deletion">\n')
    outfile.write('##ALT=<ID=INS,Description="Insertion">\n')
    outfile.write('##ALT=<ID=DUP,Description="Duplication">\n')
    outfile.write('##ALT=<ID=INV,Description="Reciprocal Inversion">\n')
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
    outfile.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV\">\n")
    outfile.write("##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strandedness\">\n")
    outfile.write("##INFO=<ID=DETAILED_TYPE,Number=1,Type=String,Description=\"Detailed type of the SV\">\n")
    outfile.write("##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Length of the unmapped sequence between breakpoint\">\n")
    outfile.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of supporting reads\">\n")
    outfile.write("##INFO=<ID=PHASESETID,Number=1,Type=String,Description=\"Matching phaseset ID for phased SVs\">\n")
    outfile.write("##INFO=<ID=HP,Number=1,Type=String,Description=\"Matching haplotype ID for phased SVs\">\n")
    outfile.write("##INFO=<ID=SUPP_READS,Number=1,Type=String,Description=\"Number of support reads\">\n")
    outfile.write("##INFO=<ID=REF_READS,Number=1,Type=String,Description=\"Number of reference reads\">\n")
    outfile.write("##INFO=<ID=CLUSTERID,Number=1,Type=String,Description=\"Cluster ID in breakpoint_graph\">\n")
    outfile.write("##INFO=<ID=INSSEQ,Number=1,Type=String,Description=\"Insertion sequence between breakpoints\">\n")
    outfile.write("##INFO=<ID=MATE_ID,Number=1,Type=String,Description=\"MATE ID for breakends\">\n")
    outfile.write("##INFO=<ID=INSIDE_VNTR,Number=1,Type=String,Description=\"True if an indel is inside a VNTR\">\n")
    outfile.write("##INFO=<ID=ALIGNED_POS,Number=1,Type=String,Description=\"Position in the reference\">\n")
    outfile.write("##INFO=<ID=BND_TYPE,Number=1,Type=String,Description=\"Junction type\">\n")
    outfile.write("##INFO=<ID=LOW_COV_IN,Number=1,Type=String,Description=\"Samples that has low coverage in that region\">\n")
    outfile.write("##INFO=<ID=INSIDE_WHITELIST,Number=1,Type=String,Description=\"SVs within the whitelist bed file\">\n")
    outfile.write("##INFO=<ID=WL_RESCUE,Number=1,Type=String,Description=\"Reason a whitelisted SV was kept despite failing a filter: WHITELIST, RECIPROCAL (a reciprocal junction corroborates it) or SUPPORTED_NEIGHBOUR (a neighbouring junction with independent read support does)\">\n")

    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    #outfile.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotyping quality\">\n")
    outfile.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of reference reads\">\n")
    outfile.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of variant reads\">\n")
    outfile.write("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n")
    outfile.write("##FORMAT=<ID=hVAF,Number=3,Type=Float,Description=\"Haplotype specific variant Allele frequency (H0,H1,H2)\">\n")
    outfile.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
    

def _sorted_breakpoints(bp_list,ref_lengths):
    db_ls = defaultdict(list)
    dbls= []
    for chr_id in ref_lengths.keys():
        db_ls[chr_id]=[]
        
    for db in bp_list:
        db_ls[db.chrom].append(db)
        
    for chr_id, dbs in db_ls.items():
        dbls += sorted(dbs, key=lambda x:x.pos)

    return dbls

    
def write_germline_vcf(vcf_list, outfile,ref_lengths):
    for db in _sorted_breakpoints(vcf_list,ref_lengths):
        outfile.write(db.to_vcf())
    outfile.close()
    
    
def write_to_vcf(double_breaks, all_ids, outpath, out_key, ref_lengths, no_ins, multisample, junction_vcf, germline_genotype, ignore_hp):            
    vcf_list = db_2_vcf(double_breaks, no_ins, all_ids, multisample, junction_vcf, germline_genotype, ignore_hp)
    key = 'somatic' if out_key == 'somatic' else 'all'
    sample_ids = [target_id.replace('.bam' , '') for target_id in all_ids]
    
    germline_outfile = open(os.path.join(outpath, 'severus_' + key + ".vcf"), "w")
    write_vcf_header(ref_lengths, germline_outfile, sample_ids)
    write_germline_vcf(vcf_list, germline_outfile,ref_lengths)
    
            
    


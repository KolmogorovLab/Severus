#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:25:41 2024

@author: keskusa2
"""

import pysam
vcf_file = 'all.delly.hg38.1kGP.ont.small.vcf.gz'
vcf = pysam.VariantFile(vcf_file)
n_samples = 1019
sv_list = []
for var in vcf:
    cc = sum([1 for idd in list(var.samples) if not var.samples[idd]['GT'] == (None, None)])/n_samples
    if var.info['SVTYPE'] == 'INS':
        sv_list.append(','.join([var.chrom, str(var.pos),var.chrom, str(var.info['SVLEN']), str(var.info['CIPOS'][1]), str(var.info['CIEND'][1]), var.info['SVTYPE'], str(cc)]))
    elif var.info['SVTYPE'] == 'BND':
        sv_list.append(','.join([var.chrom, str(var.pos),var.info['CHR2'], str(var.stop), str(var.info['CIPOS'][1]), str(var.info['CIEND'][1]), var.info['SVTYPE'], str(cc)]))
    else:
        sv_list.append(','.join([var.chrom, str(var.pos),var.chrom, str(var.stop), str(var.info['CIPOS'][1]), str(var.info['CIEND'][1]), var.info['SVTYPE'], str(cc)]))
    
out_file3 = 'PON.tsv'
with open(out_file3, "w") as fout3:
    for line in sv_list:
        fout3.write(line)
        fout3.write('\n')
fout3.close()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:05:48 2023

@author: keskusa2
"""

#!/usr/bin/env python3

import sys
import shutil
import pysam
import argparse
import os
from multiprocessing import Pool
from collections import defaultdict,Counter
import logging

from severus.build_graph import output_graphs
from severus.bam_processing import get_all_reads_parallel, init_hist, init_mm_hist, update_coverage_hist
from severus.breakpoint_finder import call_breakpoints
from severus.resolve_vntr import update_segments_by_read
from severus.__version__ import __version__


logger = logging.getLogger()


def _enable_logging(log_file, debug, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)

def _version():
    return __version__


def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.005
    MIN_BREAKPOINT_READS = 3
    MIN_MAPQ = 10
    MIN_REF_FLANK = 10000
    MAX_GENOMIC_LEN = 2000000
    MAX_SEGMENT_DIST = 1000

    #breakpoint
    BP_CLUSTER_SIZE = 50
    MIN_SV_SIZE = 50
    MIN_SV_THR = 10
    VAF_THR = 0.05
    CONTROL_VAF = 0.01
    CONTROL_COV_THR = 3
    

    SAMTOOLS_BIN = "samtools"

    parser = argparse.ArgumentParser \
        (description="Find breakpoints and build breakpoint graph from a bam file")

    parser.add_argument("-v", "--version", action="version", version=_version())
    parser.add_argument("--target-bam", dest="target_bam",
                        metavar="path", required=True, default=None, nargs="+",
                        help="path to one or multiple target bam files (e.g. tumor, must be indexed)")
    parser.add_argument("--control-bam", dest="control_bam",
                        metavar="path", required=False, default=None, nargs="+",
                        help="path to the control bam file (e.g. normal, must be indexed)")
    parser.add_argument("--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("-t", "--threads", dest="threads",
                        default=8, metavar="int", type=int, help="number of parallel threads [8]")
    parser.add_argument("--min-support", dest="bp_min_support",
                        default=0, metavar="int", type=int,
                        help=f"minimum reads supporting double breakpoint [{MIN_BREAKPOINT_READS}]")
    parser.add_argument("--min-reference-flank", dest="min_ref_flank",
                        default=MIN_REF_FLANK, metavar="int", type=int,
                        help=f"minimum distance between a breakpoint and reference ends [{MIN_REF_FLANK}]")
    parser.add_argument("--bp-cluster-size", dest="bp_cluster_size",
                        default=BP_CLUSTER_SIZE, metavar="int", type=int,
                        help=f"maximum distance in bp cluster [{BP_CLUSTER_SIZE}]")
    parser.add_argument("--min-sv-size", dest="min_sv_size",
                        default=MIN_SV_SIZE, metavar="int", type=int,
                        help=f"minimim size for sv [{MIN_SV_SIZE}]")
    parser.add_argument("--max-read-error", dest="max_read_error",
                        default=MAX_READ_ERROR, metavar="float", type=float,
                        help=f"maximum base alignment error [{MAX_READ_ERROR}]")
    parser.add_argument("--min-mapq", dest="min_mapping_quality",
                        default=MIN_MAPQ, metavar="int", type=int,
                        help=f"minimum mapping quality for aligned segment [{MIN_MAPQ}]")
    parser.add_argument("--reference-adjacencies", action="store_true", dest="reference_adjacencies",
                        default=False, help="draw reference adjacencies [False]")
    parser.add_argument("--write-alignments", action="store_true", dest="write_alignments",
                        default=False, help="write read alignments to file [False]")
    parser.add_argument("--single-bp", action="store_true", dest="single_bp",
                        default=False, help="Add hagning breakpoints [False]")
    parser.add_argument("--max-genomic-len", dest="max_genomic_len",
                        default=MAX_GENOMIC_LEN, metavar="int", type=int,
                        help=f"maximum length of genomic segment to form connected components [{MAX_GENOMIC_LEN}]")
    parser.add_argument("--phasing-vcf", dest="phase_vcf", metavar="path", help="path to vcf file used for phasing (if using haplotype specific SV calling)[None]")
    parser.add_argument("--vntr-bed", dest="vntr_file", metavar="path", help="bed file with tandem repeat locations [None]")
    parser.add_argument("--TIN-ratio", dest='control_vaf', metavar="float", type=float, default = CONTROL_VAF, help = 'Tumor in normal ratio[{CONTROL_VAF}]')
    parser.add_argument("--control-cov-thr", dest='cov_thr', metavar="float", type=int, default = CONTROL_COV_THR, help = 'Min normal coverage[{CONTROL_COV_THR}]')
    parser.add_argument("--vaf-thr", dest='vaf_thr', metavar="float", type=float, default = VAF_THR, help = 'Tumor in normal ratio[{CONTROL_VAF}]')
    parser.add_argument("--write-collapsed-dup", dest='write_segdup', action = "store_true", help = 'outputs a bed file with identified collapsed duplication regions')
    parser.add_argument("--no-ins-seq", dest='no_ins', action = "store_true", help = 'do not output insertion sequences to the vcf file')
    parser.add_argument("--output-LOH", dest='output_loh', action = "store_true", help = 'outputs a bed file with predicted LOH regions')
    parser.add_argument("--resolve-overlaps", dest='resolve_overlaps', action = "store_true")
    parser.add_argument("--ins-to-tra", dest='tra_to_ins', action = "store_false", help = 'converts insertions to translocations if mapping is known')
    parser.add_argument("--output-read-ids", dest='output_read_ids', action = "store_true", help = 'to output supporting read ids')
    parser.add_argument("--between-junction-ins", dest='ins_seq', action = "store_true", help = 'reports unmapped sequence between breakpoints')
    parser.add_argument("--max-unmapped-seq", dest='max_segment_dist',default=MAX_SEGMENT_DIST, metavar="int", type=int, help = 'maximum length of unmapped sequence between two mapped segments (if --between-junction-ins is selected the unmapped sequnce will be reported in the vcf)'')')
    parser.add_argument("--use-supplementary-tag", dest='use_supplementary_tag', action = "store_true", help = 'Uses haplotype tag in supplementary alignments')
    parser.add_argument("--PON", dest='pon_file', metavar="path", help = 'Uses PON data')
    parser.add_argument("--low-quality", dest='multisample', action = "store_true", help = 'Uses set of parameters optimized for the analysis with lower quality') 
    parser.add_argument("--target-sample", dest='target_name', default = '', help = 'Sample name for the target bams', nargs="+")
    parser.add_argument("--control-sample", dest='control_name', default = '', help = 'Sample name for the control bam', nargs="+")
    
    args = parser.parse_args()
    
    MIN_ALIGNED_LENGTH = 7000 if not args.multisample else 5000
    
    args.only_germline = False
    if args.control_bam is None:
        args.control_bam = []
        args.only_germline = True
    all_bams = args.target_bam + args.control_bam
    target_genomes = list(set(args.target_bam))
    control_genomes = list(set(args.control_bam))

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
        
    log_file = os.path.join(args.out_dir, "severus.log")
    _enable_logging(log_file, debug=False, overwrite=True)

    logger.info("Starting Severus " + _version())
    logger.debug("Cmd: %s", " ".join(sys.argv))
    logger.debug("Python version: " + sys.version)
    
    args.sv_size = max(args.min_sv_size - MIN_SV_THR, MIN_SV_THR)
    
    if not shutil.which(SAMTOOLS_BIN):
        logger.error("Error: samtools not found")
        return 1
    if not len(set(args.control_bam)) == len(args.control_bam):
        logger.error("Error: Duplicated bams are not allowed")
        return 1
    
    if not len(set(args.target_bam)) == len(args.target_bam):
        logger.error("Error: Duplicated bams are not allowed")
        return 1
        
    if len(control_genomes) > 1:
        logger.error("Error: only one control bam is allowed")
        return 1
    
    if control_genomes and control_genomes[0] in target_genomes:
        logger.error("Error: Control bam also inputted as target bam")
        return 1
        
    if args.vntr_file and not (args.vntr_file.endswith('.bed') or args.vntr_file.endswith('.bed.gz')):
        logger.error("Error: VNTR annotation file should be in bed or bed.gz format")
        return 1
    
    if args.target_name and not len(args.target_name) == len(target_genomes):
        logger.error("Error: Provide a unique sample name for each sample")
        return 1
    
    if args.bp_min_support == 0:
        args.bp_min_support = 3
    else:
        args.vaf_thr = 0

    #TODO: check that all bams have the same reference
    first_bam = all_bams[0]
    ref_lengths = None
    with pysam.AlignmentFile(first_bam, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    thread_pool = Pool(args.threads)
    
    args.write_segdups_out =''
    if args.write_segdup:
        args.write_segdups_out = open(os.path.join(args.out_dir,"severus_collaped_dup.bed"), "w")
        
    args.write_log_out = ''
    if args.output_loh:
        args.write_log_out = open(os.path.join(args.out_dir,"severus_LOH.bed"), "w")
        
    args.outpath_readqual = os.path.join(args.out_dir, "read_qual.txt")
    
    
        
    segments_by_read = []
    bam_files = defaultdict(list)
    if not args.target_name:
        args.target_name = [os.path.basename(bam_file) for bam_file in args.target_bam]
    if not args.control_name:
        args.control_name = [os.path.basename(bam_file) for bam_file in args.control_bam]
        
    genome_ids = args.target_name + args.control_name
    
    dups = [item for item, count in Counter(genome_ids).items() if count > 1]
    if dups:
        logger.info("Duplicated bam names")
        genome_ids = all_bams
        for bam_file in all_bams:
            bam_files[bam_file] = bam_file
        target_genomes = args.target_bam
        control_genomes = args.control_bam
    else:
        for i, bam_file in enumerate(all_bams):
            genome_id = genome_ids[i]
            bam_files[genome_id] = bam_file
        target_genomes = args.target_name
        control_genomes = args.control_name

    args.min_aligned_length = MIN_ALIGNED_LENGTH
    coverage_histograms = init_hist(genome_ids, ref_lengths)
    mismatch_histograms = init_mm_hist(ref_lengths)
    n90 = [MIN_ALIGNED_LENGTH]
    read_qual = defaultdict(int)
    read_qual_len = defaultdict(int)
    bg_mm = []
    for i, bam_file in enumerate(all_bams):
        genome_id = genome_ids[i]
        logger.info(f"Parsing reads from {genome_id}")
        segments_by_read_bam = get_all_reads_parallel(bam_file, thread_pool, ref_lengths, genome_id,
                                                      coverage_histograms, mismatch_histograms, n90, bg_mm,read_qual,read_qual_len,args)
        segments_by_read += segments_by_read_bam

    args.min_aligned_length = min(n90) if not args.multisample else MIN_ALIGNED_LENGTH
    logger.info('Computing read quality') 
    update_segments_by_read(segments_by_read, mismatch_histograms, bg_mm, ref_lengths,read_qual,read_qual_len, args)
    
    logger.info('Computing coverage histogram')
    update_coverage_hist(coverage_histograms,genome_ids, ref_lengths, segments_by_read, control_genomes, target_genomes, args.write_log_out)

    double_breaks = call_breakpoints(segments_by_read, ref_lengths, coverage_histograms, bam_files, genome_ids, control_genomes, thread_pool, args)
    
    output_graphs(double_breaks, coverage_histograms, thread_pool, target_genomes, control_genomes, genome_ids, ref_lengths, args)

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
from collections import defaultdict
import logging

from severus.build_graph import output_graphs
from severus.bam_processing import get_all_reads_parallel, update_coverage_hist, get_read_statistics
from severus.breakpoint_finder import call_breakpoints
from severus.resolve_vntr import update_segments_by_read



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


def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.005
    MIN_BREAKPOINT_READS = 3
    MIN_MAPQ = 10
    MIN_REF_FLANK = 10000
    MAX_GENOMIC_LEN = 50000
    MIN_ALIGNED_LENGTH = 7000

    #breakpoint
    BP_CLUSTER_SIZE = 50
    MIN_SV_SIZE = 50
    MIN_SV_THR = 10
    VAF_THR = 0.1
    CONTROL_VAF = 0.01
    

    SAMTOOLS_BIN = "samtools"

    parser = argparse.ArgumentParser \
        (description="Find breakpoints and build breakpoint graph from a bam file")

    parser.add_argument("--target-bam", dest="target_bam",
                        metavar="path", required=True, default=None, nargs="+",
                        help="path to one or multiple target bam files (must be indexed)")
    parser.add_argument("--control-bam", dest="control_bam",
                        metavar="path", required=False, default=None, nargs="+",
                        help="path to one or multiple control bam files (must be indexed)")
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
                        help=f"minimum distance between breakpoint and sequence ends [{MIN_REF_FLANK}]")
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
    parser.add_argument("--phasing-vcf", dest="phase_vcf", metavar="path", help="vcf file used for phasing [None]")
    parser.add_argument("--vntr-bed", dest="vntr_file", metavar="path", help="bed file with tandem repeat locations [None]")
    parser.add_argument("--filter-small-svs", dest='filter_small_svs', action = "store_true", help = 'filters small svs < 5000')
    parser.add_argument("--TIN-ratio", dest='control_vaf', metavar="float", type=float, default = CONTROL_VAF, help = 'Tumor in normal ratio[{CONTROL_VAF}]')
    parser.add_argument("--output-all", dest='output_all', action = "store_true")
    parser.add_argument("--write-collapsed-dup", dest='write_segdup', action = "store_true")
    parser.add_argument("--no-ins-seq", dest='no_ins', action = "store_true")
    parser.add_argument("--inbetween-ins", dest='inbetween_ins', action = "store_true")
    parser.add_argument("--only-somatic", dest='only_somatic', action = "store_true")
    parser.add_argument("--output-LOH", dest='output_loh', action = "store_true")
    parser.add_argument("--omit-resolve-overlaps", dest='resolve_overlaps', action = "store_false")
    parser.add_argument("--tra-to-ins", dest='tra_to_ins', action = "store_true")
    
    
    args = parser.parse_args()
    
    args.only_germline = False
    if args.control_bam is None:
        args.control_bam = []
        args.only_germline = True
    all_bams = args.target_bam + args.control_bam
    target_genomes = list(set(os.path.basename(b) for b in args.target_bam))
    control_genomes = list(set(os.path.basename(b) for b in args.control_bam))

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)
        
    log_file = os.path.join(args.out_dir, "severus.log")
    _enable_logging(log_file, debug=False, overwrite=True)
    logger.debug("Cmd: " + " ".join(sys.argv[1:]))
    
    args.sv_size = max(args.min_sv_size - MIN_SV_THR, MIN_SV_THR)
    
    if not shutil.which(SAMTOOLS_BIN):
        logger.error("samtools not found")
        return 1
    if args.bp_min_support == 0:
        args.bp_min_support = 3
        args.vaf_thr = VAF_THR
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
    genome_ids=[]
    n90 = [MIN_ALIGNED_LENGTH]
    for bam_file in all_bams:
        genome_id = os.path.basename(bam_file)
        genome_ids.append(genome_id)
        logger.info(f"Parsing reads from {genome_id}")
        segments_by_read_bam = get_all_reads_parallel(bam_file, thread_pool, ref_lengths, genome_id,
                                                      args.min_mapping_quality, args.sv_size)
        n90.append(get_read_statistics(segments_by_read_bam))
        segments_by_read += segments_by_read_bam
    
    args.min_aligned_length = min(n90)
    logger.info('Computing read quality') 
    update_segments_by_read(segments_by_read, ref_lengths, args)
    
    logger.info('Computing coverage histogram')
    coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, control_genomes, target_genomes, args.write_log_out)
    double_breaks = call_breakpoints(segments_by_read, ref_lengths, coverage_histograms, genome_ids, control_genomes, args)
    
    output_graphs(double_breaks, coverage_histograms, thread_pool, target_genomes, control_genomes, ref_lengths, args)

if __name__ == "__main__":
    main()
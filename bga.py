#!/usr/bin/env python3

import sys
import shutil
import pysam
import argparse
import os
from multiprocessing import Pool
from collections import defaultdict
from build_graph import build_breakpoint_graph, output_clusters_graphvis, output_clusters_csv
from breakpoint_finder import call_breakpoints, output_breaks, get_genomic_segments
from bam_processing import get_all_reads_parallel, update_coverage_hist


def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_BREAKPOINT_READS = 3
    MIN_MAPQ = 10
    #MIN_DOUBLE_BP_READS = 5
    MIN_REF_FLANK = 10000
    MAX_GENOMIC_LEN = 200000

    #breakpoint
    BP_CLUSTER_SIZE = 50
    MAX_UNALIGNED_LEN = 500
    MIN_SV_SIZE = 50
    #MIN_GRAPH_SUPPORT = 3

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
                        default=MIN_BREAKPOINT_READS, metavar="int", type=int,
                        help=f"minimum reads supporting double breakpoint [{MIN_BREAKPOINT_READS}]")
    #parser.add_argument("--double-breakpoint-min-reads", dest="double_bp_min_reads",
    #                    default=MIN_DOUBLE_BP_READS, metavar="int", type=int, help="minimum reads in double breakpoint [5]")
    parser.add_argument("--min-reference-flank", dest="min_ref_flank",
                        default=MIN_REF_FLANK, metavar="int", type=int,
                        help=f"minimum distance between breakpoint and sequence ends [{MIN_REF_FLANK}]")
    parser.add_argument("--bp_cluster_size", dest="bp_cluster_size",
                        default=BP_CLUSTER_SIZE, metavar="int", type=int,
                        help=f"maximum distance in bp cluster [{BP_CLUSTER_SIZE}]")
    parser.add_argument("--min_sv_size", dest="sv_size",
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
    parser.add_argument("--write_alignments", action="store_true", dest="write_alignments",
                        default=False, help="write read alignments to file [False]")
    parser.add_argument("--single-bp", action="store_true", dest="single_bp",
                        default=False, help="Add hagning breakpoints [False]")
    parser.add_argument("--max-genomic-len", dest="max_genomic_len",
                        default=MAX_GENOMIC_LEN, metavar="int", type=int,
                        help=f"maximum length of genomic segment to form connected components [{MAX_GENOMIC_LEN}]")
    parser.add_argument("--phasing-vcf", dest="phase_vcf", metavar="path", help="vcf file used for phasing [None]")
    parser.add_argument("--vntr-bed", dest="vntr_file", metavar="path", help="bed file with tandem repeat locations [None]")
    args = parser.parse_args()

    if args.control_bam is None:
        args.control_bam = []
    all_bams = args.target_bam + args.control_bam
    target_genomes = set(os.path.basename(b) for b in args.target_bam)
    control_genomes = set(os.path.basename(b) for b in args.control_bam)

    if not shutil.which(SAMTOOLS_BIN):
        print("samtools not found", file=sys.stderr)
        return 1

    #TODO: check that all bams have the same reference
    first_bam = all_bams[0]
    ref_lengths = None
    with pysam.AlignmentFile(first_bam, "rb") as a:
        ref_lengths = dict(zip(a.references, a.lengths))

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    thread_pool = Pool(args.threads)
    out_breakpoint_graph = os.path.join(args.out_dir, "breakpoint_graph.gv")
    out_clustered_breakpoints = os.path.join(args.out_dir, "breakpoint_clusters.csv")
    
    segments_by_read = defaultdict(list)
    genome_ids=[]
    for bam_file in all_bams:
        genome_id = os.path.basename(bam_file)
        genome_ids.append(genome_id)
        print("Parsing reads from", genome_id, file=sys.stderr)
        segments_by_read_bam =  get_all_reads_parallel(bam_file, thread_pool, ref_lengths,
                                   args.min_mapping_quality, genome_id,args.sv_size)
        segments_by_read.update(segments_by_read_bam)
        print("Parsed {0} segments".format(len(segments_by_read_bam)), file=sys.stderr)
     
    print('Computing coverage histogram')
    coverage_histograms = update_coverage_hist(genome_ids, ref_lengths, segments_by_read, args.min_mapping_quality, args.max_read_error)
    double_breaks = call_breakpoints(segments_by_read, thread_pool, ref_lengths, coverage_histograms, args)
    print('Computing segment coverage')
    genomic_segments, hb_points = get_genomic_segments(double_breaks, coverage_histograms, thread_pool, args.phase_vcf)

    print('Writing breakpoints')
    output_breaks(double_breaks, genome_ids, args.phase_vcf, open(os.path.join(args.out_dir,"breakpoints_double.csv"), "w"))
    print('Preparing graph')
    graph, adj_clusters, key_to_color = build_breakpoint_graph(double_breaks, genomic_segments, hb_points, args.max_genomic_len,
                                                                args.reference_adjacencies, target_genomes, control_genomes)
    output_clusters_graphvis(graph, adj_clusters, key_to_color, out_breakpoint_graph)
    output_clusters_csv(graph, adj_clusters, out_clustered_breakpoints)


if __name__ == "__main__":
    main()

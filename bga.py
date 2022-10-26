#!/usr/bin/env python3

import sys
import shutil
import pysam
import argparse
import os
import math
from multiprocessing import Pool

from build_graph import build_breakpoint_graph
from bam_processing import get_all_reads_parallel, get_segments_coverage
from breakpoint_finder import resolve_overlaps, get_breakpoints, get_2_breaks, DoubleBreak, output_single_breakpoints, output_breaks


def enumerate_read_breakpoints(split_reads, bp_clusters, clust_len, max_unaligned_len, bam_files, compute_coverage,
                               thread_pool, ref_lengths, out_file, min_mapq):
    """
    For each read, generate the list of breakpints
    """
    def _normalize_coord(coord, clusters):
        for cl in clusters:
            if abs(abs(coord) - cl.position) < clust_len:
                return int(math.copysign(1, coord)), cl
        return None, None

    fout = open(out_file, "w")

    #outputing adjacency connections
    for read_segments in split_reads:
        bp_names = []
        genome_id = read_segments[0].genome_id
        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):

            #TODO: consider cases with inserted sequence
            if abs(s1.read_end - s2.read_start) > max_unaligned_len or s1.mapq < min_mapq or s2.mapq < min_mapq:
                continue

            ref_bp_1 = s1.ref_end if s1.strand == "+" else s1.ref_start
            ref_bp_2 = s2.ref_start if s2.strand == "+" else s2.ref_end

            sign_1 = 1 if s1.strand == "+" else -1
            sign_2 = -1 if s2.strand == "+" else 1

            _, bp_1 = _normalize_coord(sign_1 * ref_bp_1, bp_clusters[s1.ref_id])
            _, bp_2 = _normalize_coord(sign_2 * ref_bp_2, bp_clusters[s2.ref_id])

            if not None in [bp_1, bp_2]:
                strand_1 = "+" if sign_1 > 0 else "-"
                label_1 = "{0}{1}:{2}".format(strand_1, bp_1.ref_id, bp_1.position)
                strand_2 = "+" if sign_2 > 0 else "-"
                label_2 = "{0}{1}:{2}".format(strand_2, bp_2.ref_id, bp_2.position)
                bp_names.append(label_1)
                bp_names.append(label_2)

        if len(bp_names) > 0:
            fout.write("Q {0} {1} {2}\n".format(genome_id, read_segments[0].read_id, " ".join(bp_names)))

        #output all single breakpoints
        single_bps = []
        genome_id = read_segments[0].genome_id
        for i, seg in enumerate(read_segments):
            if seg.mapq < min_mapq:
                continue

            ref_bp_1 = seg.ref_start if seg.strand == "+" else seg.ref_end
            sign_1 = -1 if seg.strand == "+" else 1
            _, bp_1 = _normalize_coord(sign_1 * ref_bp_1, bp_clusters[seg.ref_id])

            ref_bp_2 = seg.ref_end if seg.strand == "+" else seg.ref_start
            sign_2 = 1 if seg.strand == "+" else -1
            _, bp_2 = _normalize_coord(sign_2 * ref_bp_2, bp_clusters[seg.ref_id])

            if i > 0 and bp_1 is not None:
                strand_1 = "+" if sign_1 > 0 else "-"
                label_1 = "{0}{1}:{2}".format(strand_1, bp_1.ref_id, bp_1.position)
                single_bps.append(label_1)

            if i < len(read_segments) - 1 and bp_2 is not None:
                strand_2 = "+" if sign_2 > 0 else "-"
                label_2 = "{0}{1}:{2}".format(strand_2, bp_2.ref_id, bp_2.position)
                single_bps.append(label_2)

        if len(single_bps) > 0:
            fout.write("S {0} {1} {2}\n".format(genome_id, read_segments[0].read_id, " ".join(single_bps)))

    #output reference segments
    segments = []
    for seq in bp_clusters:
        clusters = bp_clusters[seq]
        if len(clusters) == 0:
            segments.append((seq, 0, ref_lengths[seq]))
        else:
            segments.append((seq, 0, clusters[0].position))

        for cl_1, cl_2 in zip(clusters[:-1], clusters[1:]):
            segments.append((seq, cl_1.position, cl_2.position))

        if len(clusters) > 0:
            segments.append((seq, clusters[-1].position, ref_lengths[seq]))

    seg_coverage = {}
    if compute_coverage:
        seg_coverage = get_segments_coverage(bam_files, segments, thread_pool)

    for seq, start, end in segments:
        label_1 = "-{0}:{1}".format(seq, start)
        label_2 = "+{0}:{1}".format(seq, end)
        coverage = 0
        if compute_coverage:
            coverage = seg_coverage[seq, start, end]

        #print(cl_1.position, cl_2.position, coverage)
        fout.write("G {0} {1} {2} {3}\n".format(label_1, label_2, end - start, int(coverage)))


def main():
    # default tunable parameters
    MAX_READ_ERROR = 0.1
    MIN_BREAKPOINT_READS = 3
    MIN_MAPQ = 10
    #MIN_DOUBLE_BP_READS = 5
    MIN_REF_FLANK = 0
    MAX_GENOMIC_LEN = 100000

    #breakpoint
    BP_CLUSTER_SIZE = 100
    MIN_SEGMENT_OVERLAP = 100
    MAX_SEGEMNT_OVERLAP = 500
    MAX_UNALIGNED_LEN = 500
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
                        default=MIN_BREAKPOINT_READS, metavar="int", type=int, help="minimum reads supporting double breakpoint [5]")
    #parser.add_argument("--double-breakpoint-min-reads", dest="double_bp_min_reads",
    #                    default=MIN_DOUBLE_BP_READS, metavar="int", type=int, help="minimum reads in double breakpoint [5]")
    #parser.add_argument("--cluster-size", dest="cluster_size",
    #                    default=BP_CLUSTER_SIZE, metavar="int", type=int, help="size of breakpoint cluster in bp [100]")
    parser.add_argument("--min-reference-flank", dest="min_ref_flank",
                        default=MIN_REF_FLANK, metavar="int", type=int,
                        help="minimum distance between breakpoint and sequence ends [0]")
    parser.add_argument("--max-read-error", dest="max_read_error",
                        default=MAX_READ_ERROR, metavar="float", type=float, help="maximum base alignment error [0.1]")
    parser.add_argument("--min-mapq", dest="min_mapping_quality",
                        default=MIN_MAPQ, metavar="int", type=int, help="minimum mapping quality for aligned segment [1]")
    #parser.add_argument("--coverage", action="store_true", dest="coverage",
    #                    default=True, help="add coverage info to breakpoint graphs")
    parser.add_argument("--reference-adjacencies", action="store_true", dest="reference_adjacencies",
                        default=False, help="draw reference adjacencies")
    parser.add_argument("--max-genomic-len", dest="max_genomic_len",
                        default=MAX_GENOMIC_LEN, metavar="int", type=int,
                        help="maximum length of genomic segment to form connected components [100000]")
    args = parser.parse_args()

    if args.control_bam is None:
        args.control_bam = []
    all_bams = args.target_bam + args.control_bam
    target_genomes = set(os.path.basename(b) for b in args.target_bam)
    control_genomes = set(os.path.basename(b) for b in args.control_bam)

    #print("target bams", args.target_bam, target_genomes, file=sys.stderr)
    #print("control bams", args.control_bam, control_genomes, file=sys.stderr)
    #print("all bams", all_bams, file=sys.stderr)

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

    out_breaks = os.path.join(args.out_dir, "breakpoints_double.csv")
    out_single_bp = os.path.join(args.out_dir, "breakpoints_single.csv")
    out_breakpoints_per_read = os.path.join(args.out_dir, "read_breakpoints")
    out_breakpoint_graph = os.path.join(args.out_dir, "breakpoint_graph.dot")

    thread_pool = Pool(args.threads)

    aln_dump_stream = open(os.path.join(args.out_dir, "read_alignments"), "w")
    all_reads = []
    split_reads = []
    for bam_file in all_bams:
        genome_tag = os.path.basename(bam_file)
        print("Parsing reads from", genome_tag, file=sys.stderr)
        genome_reads = get_all_reads_parallel(bam_file, thread_pool, aln_dump_stream,
                                              args.max_read_error, args.min_mapping_quality, genome_tag)
        all_reads.extend(genome_reads)

        genome_split_reads = [r for r in genome_reads if len(r) > 1]
        split_reads.extend(genome_split_reads)
        print("Parsed {0} reads {1} split reads".format(len(all_reads), len(split_reads)), file=sys.stderr)

    split_reads = resolve_overlaps(split_reads, MIN_SEGMENT_OVERLAP, MAX_SEGEMNT_OVERLAP)
    bp_clusters = get_breakpoints(all_reads, split_reads, BP_CLUSTER_SIZE, MAX_UNALIGNED_LEN,
                                  args.bp_min_support, args.min_ref_flank, ref_lengths, args.min_mapping_quality)
    all_breaks, balanced_breaks = get_2_breaks(bp_clusters, BP_CLUSTER_SIZE, args.bp_min_support)

    enumerate_read_breakpoints(split_reads, bp_clusters, BP_CLUSTER_SIZE, MAX_UNALIGNED_LEN, all_bams,
                               True, thread_pool, ref_lengths, out_breakpoints_per_read, args.min_mapping_quality)

    output_single_breakpoints(bp_clusters, out_single_bp)
    output_breaks(all_breaks, open(out_breaks, "w"))

    #out_balanced_breaks = os.path.join(args.out_dir, "breakpoints_balanced.csv")
    #out_inversions = os.path.join(args.out_dir, "inversions.bed")
    #output_breaks(balanced_breaks, open(out_balanced_breaks, "w"))
    #output_inversions(balanced_breaks, open(out_inversions, "w"))

    build_breakpoint_graph(out_breakpoints_per_read, args.bp_min_support, args.reference_adjacencies,
                           out_breakpoint_graph, args.max_genomic_len, target_genomes, control_genomes)


if __name__ == "__main__":
    main()

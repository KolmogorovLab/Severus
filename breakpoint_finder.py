#!/usr/bin/env python3

"""
This script takes (phased) bam file as input, and outputs coordinates
of tentative 2-breaks, in particular inversion coordinates
"""

import sys
import re
import shutil
import numpy as np
import math
from copy import copy
from collections import namedtuple, defaultdict
import pysam
from multiprocessing import Pool
import random
import os
import subprocess


class ReadConnection(object):
    __slots__ = ("ref_id_1", "pos_1", "sign_1", "ref_id_2", "pos_2", "sign_2",
                 "haplotype_1", "haplotype_2", "read_id", "genome_id")

    def __init__(self, ref_id_1, pos_1, sign_1, ref_id_2, pos_2, sign_2, haplotype_1, haplotype_2, read_id, genome_id):
        self.ref_id_1 = ref_id_1
        self.ref_id_2 = ref_id_2
        self.pos_1 = pos_1
        self.pos_2 = pos_2
        self.sign_1 = sign_1
        self.sign_2 = sign_2
        self.haplotype_1 = haplotype_1
        self.haplotype_2 = haplotype_2
        self.read_id = read_id
        self.genome_id = genome_id

    def signed_coord_1(self):
        return self.sign_1 * self.pos_1

    def signed_coord_2(self):
        return self.sign_2 * self.pos_2


class Breakpoint(object):
    __slots__ = "ref_id", "position", "spanning_reads", "connections"

    def __init__(self, ref_id, ref_position):
        self.ref_id = ref_id
        self.position = ref_position
        self.spanning_reads = []
        self.connections =[]


class DoubleBreak(object):
    __slots__ = "bp_1", "direction_1", "bp_2", "direction_2", "connections"

    def __init__(self, bp_1, direction_1, bp_2, direction_2):
        self.bp_1 = bp_1
        self.bp_2 = bp_2
        self.direction_1 = direction_1
        self.direction_2 = direction_2
        self.connections = []

    def directional_coord_1(self):
        return self.direction_1 * self.bp_1.position

    def directional_coord_2(self):
        return self.direction_2 * self.bp_2.position
        #self.spanning_reads = []

    def to_string(self):
        strand_1 = "+" if self.direction_1 > 0 else "-"
        strand_2 = "+" if self.direction_2 > 0 else "-"
        label_1 = "{0}{1}:{2}".format(strand_1, self.bp_1.ref_id, self.bp_1.position)
        label_2 = "{0}{1}:{2}".format(strand_2, self.bp_2.ref_id, self.bp_2.position)
        if label_2[1:] < label_1[1:]:
            label_1, label_2 = label_2, label_1
        bp_name = label_1 + "|" + label_2
        return bp_name


def resolve_overlaps(split_reads, min_ovlp_len, max_overlap):
    """
    Some supplementary alignments may be overlapping (e.g. in case of inversions with flanking repeat).
    This function checks if the overlap has ok structe, trims and outputs non-overlapping alignments
    """
    def _get_ovlp(seg_1, seg_2):
        max_ovlp_len = min(seg_1.read_end - seg_1.read_start, seg_2.read_end - seg_2.read_start) // 2
        if (seg_1.read_end - seg_2.read_start > min_ovlp_len and
            seg_1.read_end - seg_2.read_start < max_ovlp_len and
            seg_2.read_end > seg_1.read_end and
            seg_2.read_start > seg_1.read_start):

            return seg_1.read_end - seg_2.read_start
        else:
            return 0

    new_reads = []
    for read_segments in split_reads:
        upd_segments = []
        for i in range(len(read_segments)):
            left_ovlp = 0
            right_ovlp = 0
            if i > 0 and read_segments[i - 1].ref_id == read_segments[i].ref_id:
                left_ovlp = _get_ovlp(read_segments[i - 1], read_segments[i])
            if i < len(read_segments) - 1 and read_segments[i].ref_id == read_segments[i + 1].ref_id:
                right_ovlp = _get_ovlp(read_segments[i], read_segments[i + 1])

            left_ovlp = left_ovlp // 2
            right_ovlp = right_ovlp // 2
            seg = read_segments[i]

            if max_overlap > left_ovlp and left_ovlp > 0:
                if seg.strand == "+":
                    seg.read_start = seg.read_start + left_ovlp
                    seg.ref_start = seg.ref_start + left_ovlp
                else:
                    seg.read_start = seg.read_start + left_ovlp
                    seg.ref_end = seg.ref_end - left_ovlp
            if max_overlap > right_ovlp and right_ovlp > 0:
                if seg.strand == "+":
                    seg.read_end = seg.read_end - right_ovlp
                    seg.ref_end = seg.ref_end - right_ovlp
                else:
                    seg.read_end = seg.read_end - right_ovlp
                    seg.ref_start = seg.ref_start + right_ovlp

            upd_segments.append(seg)

        new_reads.append(upd_segments)

    return new_reads


def get_breakpoints(all_reads, split_reads, clust_len, max_unaligned_len,
                    min_reads, min_ref_flank, ref_lengths, min_mapq):
    """
    Finds regular 1-sided breakpoints, where split reads consistently connect
    two different parts of the genome
    """
    seq_breakpoints = defaultdict(list)

    def _signed_breakpoint(seg, direction):
        ref_bp, sign = None, None
        if direction == "right":
            ref_bp = seg.ref_end if seg.strand == "+" else seg.ref_start
            sign = 1 if seg.strand == "+" else -1
        elif direction == "left":
            ref_bp = seg.ref_start if seg.strand == "+" else seg.ref_end
            sign = -1 if seg.strand == "+" else 1
        return ref_bp, sign

    def _add_double(seg_1, seg_2):
        ref_bp_1, sign_1 = _signed_breakpoint(s1, "right")
        ref_bp_2, sign_2 = _signed_breakpoint(s2, "left")
        seq_breakpoints[s1.ref_id].append(ReadConnection(s1.ref_id, ref_bp_1, sign_1, s2.ref_id, ref_bp_2, sign_2,
                                                         s1.haplotype, s2.haplotype, s1.read_id, s1.genome_id))
        seq_breakpoints[s2.ref_id].append(ReadConnection(s2.ref_id, ref_bp_2, sign_2, s1.ref_id, ref_bp_1, sign_1,
                                                             s2.haplotype, s1.haplotype, s2.read_id, s2.genome_id))

    def _add_single(seg, direction):
        ref_bp, sign = _signed_breakpoint(seg, direction)
        seq_breakpoints[seg.ref_id].append(ReadConnection(seg.ref_id, ref_bp, sign, None, None, None,
                                                          seg.haplotype, None, seg.read_id, seg.genome_id))

    for read_segments in split_reads:
        #if read_segments[0].mapq >= min_mapq:
        #    _add_single(read_segments[0], "left")
        #if read_segments[-1].mapq >= min_mapq:
        #    _add_single(read_segments[-1], "right")

        for s1, s2 in zip(read_segments[:-1], read_segments[1:]):
            #TODO: consider cases with inserted sequence
            if abs(s1.read_end - s2.read_start) <= max_unaligned_len and s1.mapq >= min_mapq and s2.mapq >= min_mapq:
                _add_double(s1, s2)
            else:
                if s1.mapq >= min_mapq:
                    _add_single(s1, "right")
                if s2.mapq >= min_mapq:
                    _add_single(s2, "left")

    bp_clusters = defaultdict(list)
    for seq, bp_pos in seq_breakpoints.items():
        clusters = []
        cur_cluster = []
        bp_pos.sort(key=lambda bp: bp.pos_1)

        for rc in bp_pos:
            if cur_cluster and rc.pos_1 - cur_cluster[0].pos_1 > clust_len:
                clusters.append(cur_cluster)
                cur_cluster = [rc]
            else:
                cur_cluster.append(rc)
        if cur_cluster:
            clusters.append(cur_cluster)

        for cl in clusters:
            unique_reads = set()
            for x in cl:
                unique_reads.add((x.read_id, x.genome_id))

            by_genome_id = defaultdict(int)
            for read in unique_reads:
                by_genome_id[read[1]] += 1

            #if in at least one sample there are more than X supporting reads
            if max(by_genome_id.values()) >= min_reads:
                position = int(np.median([x.pos_1 for x in cl]))
                if position > min_ref_flank and position < ref_lengths[seq] - min_ref_flank:
                    bp_cluster = Breakpoint(seq, position)
                    bp_cluster.connections = cl
                    bp_clusters[seq].append(bp_cluster)

    #find reads that span putative breakpoints
    for read_segments in all_reads:
        for seg in read_segments:
            for bp in bp_clusters[seg.ref_id]:
                if bp.position - seg.ref_start > clust_len and seg.ref_end - bp.position > clust_len:
                    bp.spanning_reads.append(seg)

    return bp_clusters


def get_2_breaks(bp_clusters, clust_len, min_connections):
    """
    Matches one-sided breakpoints into 2-breaks, complementary pairs of one-sided breakpoints
    """
    def _normalize_coord(coord, clusters):
        for cl in clusters:
            if abs(abs(coord) - cl.position) < clust_len:
                return int(math.copysign(1, coord)), cl
        return None, None

    #Breakpoints are clustered by a single reference coordinate.
    #Our goal here is to translate it into connections (defined for two directional reference coordinates)
    #and then find balanced breakes (2-breaks)

    double_breaks = []
    double_connections = defaultdict(list)
    for seq in bp_clusters:
        for cl in bp_clusters[seq]:
            for conn in cl.connections:

                #single breakpoint
                if conn.ref_id_2 is None:
                    continue

                #normalizing the coordinates of the original read, wrt to clustered breakpoint coordinates
                dir_1, bp_1 = _normalize_coord(conn.signed_coord_1(), bp_clusters[conn.ref_id_1])
                dir_2, bp_2 = _normalize_coord(conn.signed_coord_2(), bp_clusters[conn.ref_id_2])
                if None in [bp_1, bp_2]:
                    continue

                #same breakpoint, just connected from two different sides
                if bp_1 == bp_2:
                    continue

                #new_conn = conn
                #if bp_2.position < bp_1.position:
                #    bp_1, bp_2 = bp_2, bp_1
                #    dir_1, dir_2 = dir_2, dir_1
                #    new_conn = ReadConnection(conn.ref_id_2, conn.pos_2, conn.sign_2, conn.ref_id_1,
                #                              conn.pos_1, conn.sign_1, conn.haplotype_2, conn.haplotype_1)
                #double_connections[(bp_1, dir_1, bp_2, dir_2)].append(new_conn)
                if bp_2.position >= bp_1.position:
                    double_connections[(bp_1, dir_1, bp_2, dir_2)].append(conn)

    for (bp_1, dir_1, bp_2, dir_2), conn_list in double_connections.items():
        if len(conn_list) >= min_connections:
            double_breaks.append(DoubleBreak(bp_1, dir_1, bp_2, dir_2))
            double_breaks[-1].connections = conn_list

    double_breaks.sort(key=lambda b: (b.bp_1.ref_id, b.bp_1.position, b.direction_1))

    #find breakend pairs
    all_breakends = set()
    balanced_breaks = []
    for db in double_breaks:
        all_breakends.add(db.directional_coord_1())
        all_breakends.add(db.directional_coord_2())

    for db in double_breaks:
        if -db.directional_coord_1() in all_breakends and -db.directional_coord_2() in all_breakends:
            balanced_breaks.append(db)

    return double_breaks, balanced_breaks


def output_breaks(breaks, out_stream):
    header = "#seq_id_1\tpos_1\tdirection_1\tseq_id_2\tposition_2\tdirection_2\thaplotype_1\tsupport_1\tagainst_1\thaplotype_2\tsupport_2\tagainst_2"
    out_stream.write(header + "\n")
    for br in breaks:
        by_hp_1 = defaultdict(int)
        by_hp_2 = defaultdict(int)
        #print(br.bp_1.ref_id, br.bp_1.position, br.bp_2.ref_id, br.bp_2.position)
        for conn in br.bp_1.connections:
            #print("\t", conn.ref_id_1, conn.pos_1, conn.ref_id_1, conn.pos_2, conn.haplotype_1, conn.haplotype_2)
            by_hp_1[conn.haplotype_1] += 1
        for conn in br.bp_2.connections:
            by_hp_2[conn.haplotype_1] += 1

        max_hp_1 = max(by_hp_1, key=by_hp_1.get)
        max_hp_2 = max(by_hp_2, key=by_hp_2.get)

        num_opposing_1 = 0
        for r in br.bp_1.spanning_reads:
            if r.haplotype == max_hp_1:
                num_opposing_1 += 1

        num_opposing_2 = 0
        for r in br.bp_2.spanning_reads:
            if r.haplotype == max_hp_2:
                num_opposing_2 += 1

        out_stream.write("\t".join([br.bp_1.ref_id, str(br.bp_1.position), "+" if br.direction_1 > 0 else "-",
                                   br.bp_2.ref_id, str(br.bp_2.position), "+" if br.direction_2 > 0 else "-",
                                   str(max_hp_1), str(by_hp_1[max_hp_1]), str(num_opposing_1),
                                   str(max_hp_2), str(by_hp_2[max_hp_2]), str(num_opposing_2) + "\n"]))


def output_inversions(breaks, out_stream):
    out_stream.write("#chrom\tchrmStart\tchrmEnd\thaplotype\tsupportReads\topposingReads\n")
    for br in breaks:
        if br.bp_1.ref_id == br.bp_2.ref_id and br.direction_1 > 0 and br.direction_2 > 0:
            matching = False
            for other_br in breaks:
                if (other_br.bp_1.ref_id == br.bp_1.ref_id and other_br.bp_2.ref_id == br.bp_2.ref_id and
                        other_br.bp_1.position == br.bp_1.position and other_br.bp_2.position == br.bp_2.position and
                        other_br.direction_1 < 0 and other_br.direction_2 < 0):
                    matching = True
                    break

            if matching:
                by_hp = defaultdict(int)
                for conn in br.connections:
                    by_hp[conn.haplotype_1] += 1
                max_hp = max(by_hp, key=by_hp.get)

                num_opposing = 0
                for r in br.bp_1.spanning_reads + br.bp_2.spanning_reads:
                    if r.haplotype == max_hp:
                        num_opposing += 1
                out_stream.write("\t".join([br.bp_1.ref_id, str(br.bp_1.position), str(br.bp_2.position),
                                            str(max_hp), str(by_hp[max_hp]), str(num_opposing)]) + "\n")


def output_single_breakpoints(breakpoints, filename):
    with open(filename, "w") as f:
        f.write("#chr\tposition\tHP\tsupport_reads\tspanning_reads\n")
        for seq in breakpoints:
            for bp in breakpoints[seq]:
                by_hp = defaultdict(int)
                for conn in bp.connections:
                    by_hp[conn.haplotype_1] += 1
                max_hp = max(by_hp, key=by_hp.get)

                f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(seq, bp.position, max_hp,
                                                         len(bp.connections), len(bp.spanning_reads)))

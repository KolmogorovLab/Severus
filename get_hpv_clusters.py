#!/usr/bin/env python3

import sys
from collections import namedtuple, defaultdict
from disjoint_set import DisjointSet


PafAlignment = namedtuple("PafAlignemnt", ["ref_id", "ref_start", "ref_end", "ref_len",
                                           "qry_id", "qry_start", "qry_end", "qry_len", "de_tag",
                                           "mapq", "strand"])


def parse_paf(filename):
    buffer = []
    integrations = []
    last_qry = None
    for line in open(filename, "r"):
        fields = line.strip().split()
        if fields[0] != last_qry:
            process_read(buffer, integrations)
            buffer = [fields]
            last_qry = fields[0]
        else:
            buffer.append(fields)

    return integrations


def process_read(buffer, integrations):
    alignments = []
    for fields in buffer:
        de_tag = None
        type_tag = None
        for tag in fields[12:]:
            if tag.startswith("de"):
                de_tag = tag
            if tag.startswith("tp"):
                type_tag = tag

        if type_tag[-1] == "S":
            continue

        ref_id = "HPV" if fields[5].startswith("NC") or fields[5].startswith("K") else fields[5]
        alignments.append(PafAlignment(ref_id, int(fields[7]), int(fields[8]), int(fields[6]),
                                       fields[0], int(fields[2]), int(fields[3]), int(fields[1]),
                                       de_tag, int(fields[11]), fields[4]))

    alignments = [a for a in alignments if a.mapq >= MIN_MAPQ and a.qry_end - a.qry_start > MIN_ALN]

    has_hpv = False
    has_chr = False
    left_chr = False
    right_chr = False

    alignments.sort(key = lambda a: a.qry_start)
    for a in alignments:
        if a.ref_id.startswith("HPV"):
            has_hpv = True
        if a.ref_id.startswith("chr"):
            has_chr = True

    #if has_hpv and has_chr:
    if has_hpv:
        integrations.append(alignments)

    #if has_hpv and has_chr:
    #if has_hpv:
    #    print(alignments[0].qry_id, alignments[0].qry_len, alignments[0].qry_start, alignments[-1].qry_end)
    #    print_alignments(alignments)


def print_alignments(alignments, segment_enumerator):
    for i in range(len(alignments)):
        a = alignments[i]
        gap = None
        if i < len(alignments) - 1:
            gap = alignments[i + 1].qry_start - a.qry_end
        length_diff = a.ref_end - a.ref_start - a.qry_end + a.qry_start
        segment_id = "    "
        if segment_enumerator and i > 0 and i < len(alignments) - 1:
            segment_id = segment_enumerator(a)
        #print("\t\t{0}-{1}\t{2}:{3}-{4} {5} gap:{6} diff:{7}"
        #        .format(a.qry_start, a.qry_end, a.ref_id, a.ref_start, a.ref_end, a.strand, gap, length_diff))
        print("\t\t{0}\t{1}-{2}\t{3}:{4}-{5} {6}"
                .format(segment_id, a.qry_start, a.qry_end, a.ref_id, a.ref_start, a.ref_end, a.strand))


def enumerate_inserts(alignments):
    all_starts = defaultdict(list)
    all_ends = defaultdict(list)
    for aln in alignments:
        for i in range(1, len(aln) - 1):
            all_starts[aln[i].ref_id].append(aln[i].ref_start)
            all_ends[aln[i].ref_id].append(aln[i].ref_end)

    def form_clusters(positions):
        positions.sort()
        cluster_positions = []
        cur_cluster = []
        for p in positions:
            if not cur_cluster:
                cur_cluster = [p]
            else:
                if p - cur_cluster[0] < CLUST_SIZE:
                    cur_cluster.append(p)
                else:
                    clust_center = sum(cur_cluster) // len(cur_cluster)
                    cluster_positions.append(clust_center)
                    cur_cluster = []
        if cur_cluster:
            clust_center = sum(cur_cluster) // len(cur_cluster)
            cluster_positions.append(clust_center)

        return cluster_positions

    #start_clusters = {}
    endpoint_clusters = {}
    #for seq, positions in all_starts.items():
    for seq in all_starts:
        endpoint_clusters[seq] = form_clusters(all_starts[seq] + all_ends[seq])
        #print(seq, start_clusters[seq])
    #end_clusters = {}
    #for seq, positions in all_ends.items():
    #    endpoint_clusters[seq] = form_clusters(positions)
        #print(seq, end_clusters[seq])

    def get_standard_coord(ref_id, ref_pos, is_start):
        bp_id = None

        #clusters = start_clusters if is_start else end_clusters
        for pos in endpoint_clusters[ref_id]:
            if abs(pos - ref_pos) < CLUST_SIZE:
                bp_id = pos
                break
        return bp_id

    def get_standard_coords(seg):
        """
        left_id = None
        right_id = None
        for pos in start_clusters[seg.ref_id]:
            if abs(pos - seg.ref_start) < CLUST_SIZE:
                left_id = pos
                break
        for pos in end_clusters[seg.ref_id]:
            if abs(pos - seg.ref_end) < CLUST_SIZE:
                right_id = pos
                break

        #if left_id is None:
        #    print("NONE", seg.ref_id, seg.ref_start, seg.ref_end)
        return left_id, right_id
        """
        return get_standard_coord(seg.ref_id, seg.ref_start, True), get_standard_coord(seg.ref_id, seg.ref_end, False)

    frequencies = defaultdict(int)
    for aln in alignments:
        for i in range(1, len(aln) - 1):
            frequencies[(aln[i].ref_id, *get_standard_coords(aln[i]))] += 1

    print("Frequent fragments (>=10)")
    fragment_ids = {}
    next_id = 0
    for fragment in sorted(frequencies, key=frequencies.get, reverse=True):
        fragment_ids[fragment] = next_id
        next_id += 1
        if frequencies[fragment] >= 10:
            print("\t", fragment, "\tfreq:", frequencies[fragment], "\tID:",
                  fragment[0] + "_" + str(fragment_ids[fragment]))
    print("")

    def segment_enumerator(segment):
        left, right = get_standard_coords(segment)
        return segment.ref_id + "_" + str(fragment_ids[(segment.ref_id, left, right)])
    return segment_enumerator


def print_human_hpv_clusters(integrations, segment_enumerator):
    """
    Single-linkage clustering of human-hpv integrations
    """
    def breakpoint_list(alignments):
        bp = []
        for a1, a2 in zip(alignments[:-1], alignments[1:]):
            if a1.ref_id.startswith("chr") and a2.ref_id.startswith("HPV"):
                direction = "L" if a1.strand == "+" else "R"
                bp.append((a1.ref_id, a1.ref_end if a1.strand == "+" else a1.ref_start, direction))
            if a1.ref_id.startswith("HPV") and a2.ref_id.startswith("chr"):
                direction = "R" if a2.strand == "+" else "L"
                bp.append((a2.ref_id, a2.ref_start if a2.strand == "+" else a2.ref_end, direction))
        return bp

    def common_bp(bp_list_1, bp_list_2):
        num_common = 0
        for bp1 in bp_list_1:
            for bp2 in bp_list_2:
                if bp1[0] == bp2[0] and abs(bp1[1] - bp2[1]) < 100:
                    num_common += 1
        return num_common

    ds = DisjointSet()
    for aln in integrations:
        ds.find(aln[0].qry_id)

    for aln_1 in integrations:
        bp_list_1 = breakpoint_list(aln_1)
        for aln_2 in integrations:
            bp_list_2 = breakpoint_list(aln_2)
            if common_bp(bp_list_1, bp_list_2) > 0:
                ds.union(aln_1[0].qry_id, aln_2[0].qry_id)

    clusters = defaultdict(list)
    for aln in integrations:
        if len(breakpoint_list(aln)) > 0:
            clusters[ds.find(aln[0].qry_id)].append(aln)

    for cl in clusters:
        all_breakpoints_by_chr = defaultdict(list)
        for aln in clusters[cl]:
            for bp in breakpoint_list(aln):
                all_breakpoints_by_chr[(bp[0], bp[2])].append(bp[1])

        cluster_breakpoints = []
        for (chr_id, strand), positions in all_breakpoints_by_chr.items():
            positions.sort()
            last_pos = -CLUST_SIZE
            for pos in positions:
                if pos - last_pos > CLUST_SIZE:
                    last_pos = pos
                    cluster_breakpoints.append((chr_id, pos, strand))

        print("Cluster", cluster_breakpoints)
        for aln in clusters[cl]:
            print("\t", aln[0].qry_id, aln[0].qry_len)
            print_alignments(aln, segment_enumerator)
        print("")

    print("HPV-only")
    for aln in integrations:
        if len(breakpoint_list(aln)) == 0:
            print("\t", aln[0].qry_id, aln[0].qry_len)
            print_alignments(aln, segment_enumerator)


def output_enumerated_read_segments(alignments, segment_enumerator):
    for aln in alignments:
        segments = []
        for i in range(len(aln)):
            if i > 0 and i < len(aln) - 1:
                segments.append(aln[i].strand + segment_enumerator(aln[i]))
        if len(aln) >= 3:
            print(aln[0].qry_id, " ".join(segments))


MIN_ALN = 50
MIN_MAPQ = 10
CLUST_SIZE = 50


def main():
    integrations = parse_paf(sys.argv[1])
    segment_enumerator = enumerate_inserts(integrations)
    output_enumerated_read_segments(integrations, segment_enumerator)
    #print_human_hpv_clusters(integrations, segment_enumerator)


if __name__ == "__main__":
    main()

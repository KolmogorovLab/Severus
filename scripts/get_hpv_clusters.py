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
    #if len(alignments) > 1:
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
    all_endpoints = defaultdict(list)
    for aln in alignments:
        for i in range(len(aln)):
            all_endpoints[aln[i].ref_id].append(aln[i].ref_start)
            all_endpoints[aln[i].ref_id].append(aln[i].ref_end)

        #if len(aln) >= 2:
        #    all_endpoints[aln[0].ref_id].append(aln[0].ref_end)
        #    all_endpoints[aln[-1].ref_id].append(aln[-1].ref_start)

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
                    cur_cluster = [p]

        if cur_cluster:
            clust_center = sum(cur_cluster) // len(cur_cluster)
            cluster_positions.append(clust_center)

        return cluster_positions

    endpoint_clusters = {}
    for seq in all_endpoints:
        endpoint_clusters[seq] = form_clusters(all_endpoints[seq])
        #print(seq, endpoint_clusters[seq])

    def get_standard_coord(ref_id, ref_pos):
        bp_id = None

        for pos in endpoint_clusters[ref_id]:
            if abs(pos - ref_pos) <= CLUST_SIZE:
                bp_id = pos
                break
        if bp_id is None:
            print("NONE:", ref_id, ref_pos)
        return bp_id

    def get_segment_coords(seg, end_segment=None):
        if end_segment is None:
            return get_standard_coord(seg.ref_id, seg.ref_start), get_standard_coord(seg.ref_id, seg.ref_end)

        if seg.strand == "+":
            if end_segment == "First":
                return None, get_standard_coord(seg.ref_id, seg.ref_end)
            elif end_segment == "Last":
                return get_standard_coord(seg.ref_id, seg.ref_start), None
        else:
            if end_segment == "First":
                return None, get_standard_coord(seg.ref_id, seg.ref_start)
            elif end_segment == "Last":
                return get_standard_coord(seg.ref_id, seg.ref_end), None

    frequencies = defaultdict(int)
    for aln in alignments:
        for i in range(1, len(aln) - 1):
            frequencies[(aln[i].ref_id, *get_segment_coords(aln[i]))] += 1

        if len(aln) >= 2:
            frequencies[(aln[0].ref_id, *get_segment_coords(aln[0], "First"))] += 1
            frequencies[(aln[-1].ref_id, *get_segment_coords(aln[-1], "Last"))] += 1

    print("Frequent fragments (>=10)")
    fragment_ids = {}
    next_id = 0
    for fragment in sorted(frequencies, key=frequencies.get, reverse=True):
        if fragment[1] == None and fragment[2] == None:
            continue
        if fragment[0] == "HPV" and None in fragment[1:3]:
            continue

        fragment_ids[fragment] = next_id
        next_id += 1
        if frequencies[fragment] >= 10:
            print("\t", fragment, "\tfreq:", frequencies[fragment], "\tID:",
                  fragment[0] + "_" + str(fragment_ids[fragment]))
    print("")

    def segment_enumerator(segment, end_segment=None):
        left, right = get_segment_coords(segment, end_segment)
        if not (segment.ref_id, left, right) in fragment_ids:
            return None
        return segment.ref_id + "_" + str(fragment_ids[(segment.ref_id, left, right)])

    return segment_enumerator, get_standard_coord


def print_human_hpv_clusters(integrations, segment_enumerator, min_bp_cluster):
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

        if len(clusters[cl]) >= min_bp_cluster:
            print("Cluster", len(clusters[cl]), "reads", len(cluster_breakpoints), "breakpoints")
            print("Bp:", cluster_breakpoints)
            for aln in clusters[cl]:
                print("\t", aln[0].qry_id, aln[0].qry_len)
                print_alignments(aln, segment_enumerator)
            print("")

    #print("HPV-only")
    #for aln in integrations:
    #    if len(breakpoint_list(aln)) == 0:
    #        print("\t", aln[0].qry_id, aln[0].qry_len)
    #        print_alignments(aln, segment_enumerator)


def output_enumerated_read_segments(alignments, segment_enumerator):
    for aln in alignments:
        segments = []

        if len(aln) >= 2:
            first_id = segment_enumerator(aln[0], "First")
            if first_id is not None:
                segments.append(aln[0].strand + first_id)

        for i in range(1, len(aln) - 1):
            segments.append(aln[i].strand + segment_enumerator(aln[i]))

        if len(aln) >= 2:
            last_id = segment_enumerator(aln[-1], "Last")
            if last_id is not None:
                segments.append(aln[-1].strand + last_id)

        if len(segments) >= 2:
            print(aln[0].qry_id, " ".join(segments))


def _neg_sign(sign):
    if sign == "+":
        return "-"
    return "+"


def enumerate_breakpoints(alignments, coordinate_enumerator):
    for aln in alignments:
        if len(aln) == 1:
            continue
        #segments = []
        #for seg in aln:
        #    segments.append("{0}{1}_{2}:{3}".format(seg.strand, seg.ref_id, seg.ref_start, seg.ref_end))
        #print(aln[0].qry_id, " ".join(segments))

        breakpoints = []
        for s1, s2 in zip(aln[:-1], aln[1:]):
            #print(s1, s2)
            ref_bp_1 = s1.ref_end if s1.strand == "+" else s1.ref_start
            ref_bp_2 = s2.ref_start if s2.strand == "+" else s2.ref_end

            coord_1 = coordinate_enumerator(s1.ref_id, ref_bp_1)
            coord_2 = coordinate_enumerator(s2.ref_id, ref_bp_2)

            if coord_1 is not None and coord_2 is not None:
                label_1 = "{0}{1}_{2}".format(s1.strand, s1.ref_id, coord_1)
                label_2 = "{0}{1}_{2}".format(_neg_sign(s2.strand), s2.ref_id, coord_2)
                if label_2[1:] < label_1[1:]:
                    label_1, label_2 = label_2, label_1
                bp_name = label_1 + "|" + label_2
                breakpoints.append(bp_name)
            #else:
            #    print("AAA")

        if len(breakpoints) >= 1:
            print(aln[0].qry_id, " ".join(breakpoints))
        #print("")

MIN_ALN = 50
MIN_MAPQ = 10
CLUST_SIZE = 50
MIN_BP_CLUSTER = 2

def main():
    integrations = parse_paf(sys.argv[1])
    segment_enumerator, coordinate_enumerator = enumerate_inserts(integrations)
    enumerate_breakpoints(integrations, coordinate_enumerator)
    #output_enumerated_read_segments(integrations, segment_enumerator)
    #print_human_hpv_clusters(integrations, segment_enumerator, MIN_BP_CLUSTER)


if __name__ == "__main__":
    main()

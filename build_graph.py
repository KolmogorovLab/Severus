#!/usr/bin/env python3

import argparse
import sys
import networkx as nx


def neg_segment(segment):
    if segment.startswith("+"):
        return "-" + segment[1:]
    else:
        return "+" + segment[1:]


#COLORS = ["yellowgreen", "thistle", "peachpuff", "yellow", "khaki", "steelblue", "hotpink", "preu"]
COLORS = ["red", "green", "blue", "orange", "purple", "cyan", "hotpink", "preu"]


def build_graph(read_segments, kmer_size, min_coverage, max_genomic_len, reference_adjacencies):
    SEQUENCE_KEY = "__genomic"
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]

    def update_graph(segments, genome_id):
        for i in range(0, len(segments) - kmer_size):
            left_kmer = node_to_id(",".join(segments[i : i + kmer_size]))
            right_kmer = node_to_id(",".join(segments[i + 1 : i + kmer_size + 1]))

            #adjacency edges are dashed
            edge_style = "solid"
            if i % 2 == 0:
                edge_style = "dashed"

            if not g.has_node(left_kmer):
                g.add_node(left_kmer)
            if not g.has_node(right_kmer):
                g.add_node(right_kmer)

            if g.has_edge(left_kmer, right_kmer, key=genome_id):
                g[left_kmer][right_kmer][genome_id]["weight"] = g[left_kmer][right_kmer][genome_id]["weight"] + 1
            else:
                g.add_edge(left_kmer, right_kmer, key=genome_id, weight=1, style=edge_style)

    genome_segments = []
    for line in open(read_segments, "r"):
        if line.startswith("Q"):
            fields = line.strip().split()
            genome_id, read_id, segments = fields[1], fields[2], fields[3:]
            if len(segments) > kmer_size:
                update_graph(segments, genome_id)
                #reversed_segments = []
                #for seg in segments[::-1]:
                #    reversed_segments.append(neg_segment(seg))
                #update_graph(reversed_segments)

        if line.startswith("G"):
            fields = line.strip().split()
            genome_segments.append(fields[1:])

    #add sequence segments
    for gs in genome_segments:
        if int(gs[2]) < max_genomic_len:
            label="\"L:{0} C:{1}\"".format(gs[2], gs[3])
            g.add_edge(node_to_id(gs[0]), node_to_id(gs[1]), label=label, key=SEQUENCE_KEY)

    #remove edges with low coverage
    edges_to_delete = set()
    for u, v, _key in g.edges:
       if _key != SEQUENCE_KEY and g[u][v][_key]["weight"] < min_coverage:
           edges_to_delete.add((u, v, _key))
    for u, v, _key in edges_to_delete:
        g.remove_edge(u, v, key=_key)

    #remove isolated nodes
    g.remove_nodes_from(list(nx.isolates(g)))

    #add reference adjacencies if requested
    if reference_adjacencies:
        for n in node_ids:
            if n.startswith("+"):
                other_n = "-" + n[1:]
                if other_n in node_ids and g.has_node(node_ids[n]) and g.has_node(node_ids[other_n]):
                    g.add_edge(node_ids[n], node_ids[other_n], style="dashed", key=SEQUENCE_KEY)

    #adding labels to query edges
    for u, v, _key in g.edges:
        if _key != SEQUENCE_KEY:
            g[u][v][_key]["label"] = "\"R:{0}\"".format(g[u][v][_key]["weight"])

    #adding colors
    key_to_color = {}
    next_color = 0
    for u, v, _key in g.edges:
        if _key != SEQUENCE_KEY:
            if _key not in key_to_color:
                key_to_color[_key] = COLORS[next_color]
                next_color = (next_color + 1) % len(COLORS)
            g[u][v][_key]["color"] = key_to_color[_key]
    print(key_to_color)

    #label nodes
    for n in g.nodes:
        g.nodes[n]["label"] = "\\n".join(id_to_kmers[n].split(","))

    #legend edges
    for i, (key, color) in enumerate(key_to_color.items()):
        node_1 = "legend_{0}_1".format(i)
        node_2 = "legend_{0}_2".format(i)
        g.add_node(node_1, label="")
        g.add_node(node_2, label="")
        g.add_edge(node_1, node_2, color=color, label=key)

    return g


def build_breakpoint_graph(reads_segments_path, min_support, reference_adjacencies, out_file):
    KMER = 1
    MAX_GENOMIC_LEN = 1000000000

    graph = build_graph(reads_segments_path, KMER, min_support, MAX_GENOMIC_LEN, reference_adjacencies)
    nx.drawing.nx_pydot.write_dot(graph, out_file)


"""
def main():
    parser = argparse.ArgumentParser \
        (description="Build breakpoint graph")

    parser.add_argument("--reads", dest="reads_path",
                        metavar="path", required=True,
                        help="path to read breakpints file")
    parser.add_argument("--out", dest="out_graph",
                        default=None, required=True,
                        metavar="path", help="Output graph")
    parser.add_argument("--min-coverage", dest="min_coverage",
                        metavar="int", type=int, required=True, default=2,
                        help="Minimum read coverage for breakpoint edge")
    args = parser.parse_args()

    graph = build_graph(args.reads_path, KMER, args.min_coverage, MAX_GENOMIC_LEN)
    nx.drawing.nx_pydot.write_dot(graph, args.out_graph)


if __name__ == "__main__":
    main()
"""

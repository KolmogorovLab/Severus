#!/usr/bin/env python3

import sys
import networkx as nx


def neg_segment(segment):
    if segment.startswith("+"):
        return "-" + segment[1:]
    else:
        return "+" + segment[1:]

def build_graph(read_segments, kmer_size, min_coverage, max_genomic_len):
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]

    def update_graph(segments):
        for i in range(0, len(segments) - kmer_size):
            left_kmer = node_to_id(",".join(segments[i : i + kmer_size]))
            right_kmer = node_to_id(",".join(segments[i + 1 : i + kmer_size + 1]))

            if not g.has_node(left_kmer):
                g.add_node(left_kmer)
            if not g.has_node(right_kmer):
                g.add_node(right_kmer)

            if g.has_edge(left_kmer, right_kmer):
                g[left_kmer][right_kmer][0]["weight"] = g[left_kmer][right_kmer][0]["weight"] + 1
                #print(g[left_kmer][right_kmer]["weight"])
            else:
                g.add_edge(left_kmer, right_kmer, weight=1, color="red")
                #print("new node", left_kmer, right_kmer)

            #print(left_kmer, right_kmer)

    genome_segments = []
    for line in open(read_segments, "r"):
        if line.startswith("Q"):
            fields = line.strip().split()
            segments = fields[2:]
            if len(segments) <= kmer_size:
                continue

            update_graph(segments)
            #reversed_segments = []
            #for seg in segments[::-1]:
            #    reversed_segments.append(neg_segment(seg))
            #update_graph(reversed_segments)

        if line.startswith("G"):
            fields = line.strip().split()
            genome_segments.append(fields[1:])

    edges_to_delete = set()
    for u, v, _num in g.edges:
       g[u][v][_num]["label"] = str(g[u][v][_num]["weight"])
       if g[u][v][_num]["weight"] < min_coverage:
           edges_to_delete.add((u, v))
    for u, v in edges_to_delete:
        g.remove_edge(u, v)

    for gs in genome_segments:
        if int(gs[2]) < max_genomic_len:
            g.add_edge(node_to_id(gs[0]), node_to_id(gs[1]), label=gs[2])

    for n in g.nodes:
        g.nodes[n]["label"] = "\\n".join(id_to_kmers[n].split(","))

    g.remove_nodes_from(list(nx.isolates(g)))

    return g


def main():
    KMER = 1
    MIN_COVERAGE = 2
    MAX_GENOMIC_LEN = 1000000000

    graph = build_graph(sys.argv[1], KMER, MIN_COVERAGE, MAX_GENOMIC_LEN)
    nx.drawing.nx_pydot.write_dot(graph, sys.argv[2])


if __name__ == "__main__":
    main()

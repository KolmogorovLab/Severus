#!/usr/bin/env python3

import sys
import networkx as nx


def neg_segment(segment):
    if segment.startswith("+"):
        return "-" + segment[1:]
    else:
        return "+" + segment[1:]

def build_graph(read_segments, kmer_size, min_coverage):
    g = nx.DiGraph()

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
                g[left_kmer][right_kmer]["weight"] = g[left_kmer][right_kmer]["weight"] + 1
                #print(g[left_kmer][right_kmer]["weight"])
            else:
                g.add_edge(left_kmer, right_kmer, weight=1)
                #print("new node", left_kmer, right_kmer)

            #print(left_kmer, right_kmer)

    for line in open(read_segments, "r"):
        fields = line.split()
        segments = fields[1:]
        if len(segments) <= kmer_size:
            continue

        reversed_segments = []
        for seg in segments[::-1]:
            reversed_segments.append(neg_segment(seg))

        update_graph(segments)
        update_graph(reversed_segments)

    edges_to_delete = set()
    for n in g.nodes:
        g.nodes[n]["label"] = "\\n".join(id_to_kmers[n].split(","))

    for u, v in g.edges:
       g[u][v]["label"] = str(g[u][v]["weight"])
       if g[u][v]["weight"] < min_coverage:
           edges_to_delete.add((u, v))

    for u, v in edges_to_delete:
        g.remove_edge(u, v)
    g.remove_nodes_from(list(nx.isolates(g)))

    return g


def main():
    KMER = 1
    MIN_COVERAGE = 5
    graph = build_graph(sys.argv[1], KMER, MIN_COVERAGE)
    nx.drawing.nx_pydot.write_dot(graph, sys.argv[2])


if __name__ == "__main__":
    main()

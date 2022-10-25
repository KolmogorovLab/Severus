#!/usr/bin/env python3

import argparse
import sys
import os
import networkx as nx
from collections import defaultdict


def neg_segment(segment):
    if segment.startswith("+"):
        return "-" + segment[1:]
    else:
        return "+" + segment[1:]


#COLORS = ["yellowgreen", "thistle", "peachpuff", "yellow", "khaki", "steelblue", "hotpink", "preu"]
COLORS = ["red", "green", "blue", "orange", "purple", "cyan", "hotpink", "preu"]


def build_graph(read_segments, min_coverage, max_genomic_len, reference_adjacencies):
    kmer_size = 1
    SEQUENCE_KEY = "__genomic"
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]

    def update_adjacencies(segments, genome_id):
        for i in range(0, len(segments) - kmer_size):
            left_kmer = node_to_id(",".join(segments[i : i + kmer_size]))
            right_kmer = node_to_id(",".join(segments[i + 1 : i + kmer_size + 1]))

            #adjacency edges are dashed
            adjacency = (i % 2 == 0)
            if not adjacency:
                continue
            edge_style = "dashed"

            if not g.has_node(left_kmer):
                g.add_node(left_kmer)
            if not g.has_node(right_kmer):
                g.add_node(right_kmer)

            if g.has_edge(left_kmer, right_kmer, key=genome_id):
                g[left_kmer][right_kmer][genome_id]["support"] = g[left_kmer][right_kmer][genome_id]["support"] + 1
            else:
                g.add_edge(left_kmer, right_kmer, key=genome_id, support=1, style=edge_style)

    genome_segments = []
    breakpoints_coverage = defaultdict(int)

    for line in open(read_segments, "r"):
        if line.startswith("Q"):
            fields = line.strip().split()
            genome_id, read_id, segments = fields[1], fields[2], fields[3:]
            if len(segments) > kmer_size:
                update_adjacencies(segments, genome_id)

        if line.startswith("S"):
            fields = line.strip().split()
            genome_id, read_id, breakpoints = fields[1], fields[2], fields[3:]
            for bp in breakpoints:
                breakpoints_coverage[(bp, genome_id)] += 1

        if line.startswith("G"):
            fields = line.strip().split()
            genome_segments.append(fields[1:])

    #add sequence segments
    for i, gs in enumerate(genome_segments):
        if not max_genomic_len or int(gs[2]) < max_genomic_len:
            label="L:{0} C:{1}".format(gs[2], gs[3])
            g.add_edge(node_to_id(gs[0]), node_to_id(gs[1]), label=label, key=SEQUENCE_KEY)

        #break connected components
        else:
            node_1 = "broken_{0}_1".format(i)
            node_2 = "broken_{0}_2".format(i)
            g.add_node(node_1, label=gs[1], style="filled", fillcolor="grey")
            g.add_node(node_2, label=gs[0], style="filled", fillcolor="grey")

            label="L:{0} C:{1}".format(gs[2], gs[3])
            g.add_edge(node_to_id(gs[0]), node_1, label=label, key=SEQUENCE_KEY)
            g.add_edge(node_2, node_to_id(gs[1]), label=label, key=SEQUENCE_KEY)

    #add hanging segments (single breakpoints), if there is no double breakpoint
    hanging_num = 0
    for (node, genome_id), multiplicity in breakpoints_coverage.items():
        if multiplicity < min_coverage:
            continue

        has_adjacency = False
        for u, v, _key in g.edges(node_to_id(node), keys=True):
            if _key == genome_id and g[u][v][_key]["support"] >= min_coverage:
                has_adjacency = True

        if not has_adjacency:
            new_node = "hanging_{0}".format(hanging_num)
            hanging_num += 1
            g.add_node(new_node, label="", shape="point")
            g.add_edge(node_to_id(node), new_node, key=genome_id, support=multiplicity, style="dashed")

    #remove edges with low coverage
    edges_to_delete = set()
    safe_edges = set()
    for u, v, _key in g.edges:
        if _key == SEQUENCE_KEY:
            continue

        if g[u][v][_key]["support"] < min_coverage:
            edges_to_delete.add((u, v, _key))
        else:
            safe_edges.add((u, v))
            safe_edges.add((v, u))

    for u, v, _key in edges_to_delete:
        if (u, v) not in safe_edges:
            g.remove_edge(u, v, key=_key)

    #remove isolated nodes
    #g.remove_nodes_from(list(nx.isolates(g)))

    #add reference adjacencies if requested
    ref_style = "dashed" if reference_adjacencies else "invis"
    for n in node_ids:
        if n.startswith("+"):
            other_n = "-" + n[1:]
            if other_n in node_ids and g.has_node(node_ids[n]) and g.has_node(node_ids[other_n]):
                g.add_edge(node_ids[n], node_ids[other_n], style=ref_style, weight="0",
                           contraint="false", key=SEQUENCE_KEY)

    #adding labels to query edges
    for u, v, _key in g.edges:
        if _key != SEQUENCE_KEY:
            g[u][v][_key]["label"] = "R:{0}".format(g[u][v][_key]["support"])

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
        if n in id_to_kmers:
            g.nodes[n]["label"] = "\\n".join(id_to_kmers[n].split(","))

    return g, key_to_color


def add_legend(key_to_color, fout):
    fout.write("subgraph cluster_01 {\n\tlabel = \"Legend\";\n\tnode [shape=point]\n{\n\trank=same\n")

    for i, (key, color) in enumerate(key_to_color.items()):
        node_1 = "legend_{0}_1".format(i)
        node_2 = "legend_{0}_2".format(i)
        fout.write("\t{0} -- {1} [color={2}, label=\"{3}\"];\n".format(node_1, node_2, color, key))
    fout.write("}\n};\n")


def output_connected_components(graph, out_stream, target_genomes, control_genomes):
    #compute target-normal rank
    components_list = []
    for cc in nx.connected_components(graph):
        target_adj = set()
        control_adj = set()
        for u, v, key in graph.edges(cc, keys=True):
            if str(u).startswith("hanging") or str(v).startswith("hanging"):
                continue

            if key in target_genomes:
                target_adj.add((u, v))
            if key in control_genomes:
                control_adj.add((u, v))

        target_adj = len(target_adj)
        control_adj = len(control_adj)
        rank = None
        if target_adj > 0 and control_adj == 0:
            rank = target_adj + 100
        elif target_adj == 0 and control_adj == 0:
            rank = -100
        elif target_adj == 0 and control_adj > 0:
            rank = -control_adj
        else:
            rank = target_adj / control_adj

        components_list.append((rank, cc, target_adj, control_adj))

    #output in rank order as subclusters
    components_list.sort(key=lambda p: p[0], reverse=True)
    for subgr_num, (rank, cc, target_adj, control_adj) in enumerate(components_list):
        out_stream.write("subgraph cluster{0} {{\n".format(subgr_num))
        #print(subgr_num, target_adj, control_adj, rank)

        for n in cc:
            properties = []
            for key in graph.nodes[n]:
                properties.append("{0}=\"{1}\"".format(key, graph.nodes[n][key]))
            out_stream.write("{0} [{1}];\n".format(n, ",".join(properties)))

        for u, v, key in graph.edges(cc, keys=True):
            properties = []
            for prop in graph[u][v][key]:
                properties.append("{0}=\"{1}\"".format(prop, graph[u][v][key][prop]))
            out_stream.write("{0} -- {1} [{2}];\n".format(u, v, ",".join(properties)))

        out_stream.write("}\n\n")


def build_breakpoint_graph(reads_segments_path, min_support, reference_adjacencies,
                           out_file, max_genomic_len, target_genomes, control_genomes):
    graph, key_to_color = build_graph(reads_segments_path, min_support, max_genomic_len, reference_adjacencies)

    with open(out_file, "w") as out_stream:
        out_stream.write("graph {\n")
        out_stream.write("edge [penwidth=2];\n")

        add_legend(key_to_color, out_stream)
        output_connected_components(graph, out_stream, target_genomes, control_genomes)

        out_stream.write("}\n")

    #tmp_graph = out_file + "_tmp"
    #nx.drawing.nx_pydot.write_dot(graph, tmp_graph)
    #os.remove(tmp_graph)

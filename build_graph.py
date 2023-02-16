#!/usr/bin/env python3

import argparse
import sys
import os
import networkx as nx
from collections import defaultdict


#def neg_segment(segment):
#    if segment.startswith("+"):
#        return "-" + segment[1:]
#    else:
#        return "+" + segment[1:]

#COLORS = ["yellowgreen", "thistle", "peachpuff", "yellow", "khaki", "steelblue", "hotpink", "preu"]
COLORS = ["red", "green", "blue", "orange", "purple", "cyan", "hotpink", "preu"]
SEQUENCE_KEY = "__genomic"


def build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies):
    kmer_size = 1
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]

    def conv_dir(dir1):
        return '-' if dir1==-1 else '+'

    def update_adjacencies(double_bp):
        left_kmer = node_to_id(":".join([conv_dir(double_bp.direction_1)+double_bp.bp_1.ref_id, str(double_bp.bp_1.position)]))
        right_kmer = node_to_id(":".join([conv_dir(double_bp.direction_2)+double_bp.bp_2.ref_id, str(double_bp.bp_2.position)]))
        if not g.has_node(left_kmer):
            g.add_node(left_kmer,label = id_to_kmers[left_kmer])
        if not g.has_node(right_kmer):
            g.add_node(right_kmer, label = id_to_kmers[right_kmer])
        if g.has_edge(left_kmer, right_kmer, key=double_bp.genome_id):
            g[left_kmer][right_kmer][double_bp.genome_id]["support"] += double_bp.supp
        else:
            edge_weight = 2 if double_bp.genotype == "hom" else 1
            g.add_edge(left_kmer, right_kmer, key=double_bp.genome_id, support=double_bp.supp,
                       style=double_bp.edgestyle, penwidth=edge_weight, adjacency=True, genotype=double_bp.genotype)
        g[left_kmer][right_kmer][double_bp.genome_id]["label"] = "R:{0}".format(g[left_kmer][right_kmer][double_bp.genome_id]["support"])

    def add_genomic_edge(seg, ref_style):
        edge_flag = True
        left_kmer = node_to_id(":".join([seg.dir1+seg.ref_id, str(seg.pos1)]))
        right_kmer = node_to_id(":".join([seg.dir2+seg.ref_id, str(seg.pos2)]))
        left_kmer_2 = node_to_id(":".join([seg.dir2+seg.ref_id, str(seg.pos1)]))
        right_kmer_2 = node_to_id(":".join([seg.dir1+seg.ref_id, str(seg.pos2)]))

        if not g.has_node(left_kmer_2):
            g.add_node(left_kmer_2, label = id_to_kmers[left_kmer_2])
        if not g.has_node(left_kmer):
            g.add_node(left_kmer, label = id_to_kmers[left_kmer])
        g.add_edge(left_kmer, left_kmer_2, key=SEQUENCE_KEY, style=ref_style, adjacency=False)
        if not g.has_node(right_kmer_2):
            g.add_node(right_kmer_2, label = id_to_kmers[right_kmer_2])
        if not g.has_node(right_kmer):
            g.add_node(right_kmer, label = id_to_kmers[right_kmer])
        g.add_edge(right_kmer, right_kmer_2, key=SEQUENCE_KEY, style=ref_style, adjacency=False)

        if seg.pos1 in hb_points[seg.ref_id]:
            hb_attr = {left_kmer: {'style':"filled", 'fillcolor':"grey", 'label' : 'HB'},
                       left_kmer_2:{'style':"filled", 'fillcolor':"grey", 'label':'HB'}}
            nx.set_node_attributes(g, hb_attr )
        if seg.pos2 in hb_points[seg.ref_id]:
            hb_attr = {right_kmer:{'style':"filled", 'fillcolor':"grey", 'label': 'HB'},
                       right_kmer_2:{'style':"filled", 'fillcolor':"grey", 'label':'HB'}}
            nx.set_node_attributes(g, hb_attr )

        if (g.has_edge(left_kmer, right_kmer, key = seg.genome_id) and g[left_kmer][right_kmer][seg.genome_id]["style"] == 'dashed'):
            edge_flag = False
        elif g.has_edge(left_kmer_2, right_kmer, key = seg.genome_id) and g[left_kmer_2][right_kmer][seg.genome_id]["style"] == 'dashed':
            edge_flag = False
        elif g.has_edge(left_kmer, right_kmer_2, key = seg.genome_id) and g[left_kmer][right_kmer_2][seg.genome_id]["style"] == 'dashed':
            edge_flag = False
        elif g.has_edge(left_kmer_2, right_kmer_2, key = seg.genome_id)and g[left_kmer_2][right_kmer_2][seg.genome_id]["style"] == 'dashed':
            edge_flag = False

        if edge_flag:
            if g.has_edge(left_kmer, right_kmer, key = seg.genome_id):
                g[left_kmer][right_kmer][seg.genome_id]["support"] += int(seg.coverage[0])
                g[left_kmer][right_kmer][seg.genome_id]["penwidth"] = 2
            else:
                g.add_edge(left_kmer, right_kmer, key=seg.genome_id,
                           support=int(seg.coverage[0]), style='solid', penwidth = 1, adjacency=False)

            label="L:{0} C:{1}".format(seg.length_bp, g[left_kmer][right_kmer][seg.genome_id]["support"])
            g[left_kmer][right_kmer][seg.genome_id]["label"] = label

    for double_bp in double_breaks:
        update_adjacencies(double_bp)

    ref_style = "dashed" if reference_adjacencies else "invis"
    for seg in genomicsegments:
        if seg.length_bp < max_genomic_len:
            if not (seg.haplotype ==0 and seg.coverage ==0):
                add_genomic_edge(seg, ref_style)

    #remove edges with low coverage

    #adding colors
    key_to_color = {}
    next_color = 0
    for u, v, _key in g.edges:
        if _key != SEQUENCE_KEY:
            if _key not in key_to_color:
                key_to_color[_key] = COLORS[next_color]
                next_color = (next_color + 1) % len(COLORS)
            g[u][v][_key]["color"] = key_to_color[_key]
    #print(key_to_color)

    return g, key_to_color


def output_clusters_graphvis(graph, connected_components, key_to_color, out_file):
    def _add_legend(key_to_color, fout):
        fout.write("subgraph cluster_01 {\n\tlabel = \"Legend\";\n\tnode [shape=point]\n{\n\trank=same\n")
        for i, (key, color) in enumerate(key_to_color.items()):
            node_1 = "legend_{0}_1".format(i)
            node_2 = "legend_{0}_2".format(i)
            fout.write("\t{0} -- {1} [color={2}, label=\"{3}\"];\n".format(node_1, node_2, color, key))
        fout.write("}\n};\n")

    def _draw_components(components_list, out_stream):
        for subgr_num, (rank, cc, target_adj, control_adj) in enumerate(components_list):
            out_stream.write("subgraph cluster{0} {{\n".format(subgr_num))
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

    with open(out_file, "w") as out_stream:
        out_stream.write("graph {\n")
        out_stream.write("edge [penwidth=2];\n")
        _add_legend(key_to_color, out_stream)
        _draw_components(connected_components, out_stream)
        out_stream.write("}\n")


def output_clusters_csv(graph, connected_components, out_file):
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\tadj_1\tadj_2\tgenome_ids\tread_support\tgenotype\n")
        for subgr_num, (rank, cc, target_adj, control_adj) in enumerate(connected_components):
            all_adj = defaultdict(list)
            for u, v, key in graph.edges(cc, keys=True):
                if key != SEQUENCE_KEY and graph[u][v][key]["adjacency"] == True:
                    all_adj[(u, v)].append(key)

            for (u, v), keys in all_adj.items():
                label_1 = graph.nodes[u]["label"]
                label_2 = graph.nodes[v]["label"]
                keys_text = ",".join(keys)
                read_support = ",".join([str(graph[u][v][k]["support"]) for k in keys])
                genotypes = ",".join([graph[u][v][k]["genotype"] for k in keys])
                fout.write(f"bga_{subgr_num}\t{label_1}\t{label_2}\t{keys_text}\t{read_support}\t{genotypes}\n")


def cluster_adjacencies(graph, target_genomes, control_genomes):
    components_list = []
    for cc in nx.connected_components(graph):
        target_adj = set()
        control_adj = set()
        for u, v, key in graph.edges(cc, keys=True):
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
    components_list.sort(key=lambda p: p[0], reverse=True)

    return components_list


def build_breakpoint_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies,
                           target_genomes, control_genomes):
    graph, key_to_color = build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies)
    adj_clusters = cluster_adjacencies(graph, target_genomes, control_genomes)
    return graph, adj_clusters, key_to_color

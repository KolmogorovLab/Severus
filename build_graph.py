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


def build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies):
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
            g[left_kmer][right_kmer][double_bp.genome_id]["support"] = g[left_kmer][right_kmer][double_bp.genome_id]["support"] + double_bp.supp
        else:
            g.add_edge(left_kmer, right_kmer, key=double_bp.genome_id, support=double_bp.supp, style=double_bp.edgestyle, penwidth = double_bp.edgeweight)
        g[left_kmer][right_kmer][double_bp.genome_id]["label"] = "R:{0}".format(g[left_kmer][right_kmer][double_bp.genome_id]["support"])
        

    def add_genomic_edge(seg, ref_style):
        edge_flag = True
        left_kmer = node_to_id(":".join([seg.dir1+seg.ref_id, str(seg.pos1)]))
        right_kmer = node_to_id(":".join([seg.dir2+seg.ref_id, str(seg.pos2)]))
        left_kmer_2 = node_to_id(":".join([seg.dir2+seg.ref_id, str(seg.pos1)]))
        right_kmer_2 = node_to_id(":".join([seg.dir1+seg.ref_id, str(seg.pos2)]))
        ln = [ln for ln in [left_kmer, left_kmer_2] if g.has_node(ln)]
        rn=[rn for rn in [right_kmer, right_kmer_2] if g.has_node(rn)]
        if g.has_node(left_kmer) and not g.has_node(left_kmer_2):
            g.add_node(left_kmer_2, label = id_to_kmers[left_kmer_2])
            g.add_edge(left_kmer, left_kmer_2, key=SEQUENCE_KEY, style=ref_style)
        elif g.has_node(left_kmer_2) and not g.has_node(left_kmer):
            if seg.pos1 in hb_points[seg.ref_id]:
                g.add_node(left_kmer,style="filled", fillcolor="grey", label = 'HB')
            else:
                g.add_node(left_kmer, label = id_to_kmers[left_kmer])
            g.add_edge(left_kmer, left_kmer_2, key=SEQUENCE_KEY, style=ref_style)
        if g.has_node(right_kmer) and not g.has_node(right_kmer_2):
            g.add_node(right_kmer_2, label = id_to_kmers[right_kmer_2])
            g.add_edge(right_kmer, right_kmer_2, key=SEQUENCE_KEY, style=ref_style)
        elif g.has_node(right_kmer_2) and not g.has_node(right_kmer):
            if seg.pos2 in hb_points[seg.ref_id]:
                g.add_node(right_kmer,style="filled", fillcolor="grey", label = 'HB')
            else:
                g.add_node(right_kmer, label = id_to_kmers[right_kmer])
            g.add_edge(right_kmer, right_kmer_2, key=SEQUENCE_KEY, style=ref_style)
        if (g.has_edge(left_kmer, right_kmer, key = seg.genome_id) and g[left_kmer][right_kmer][seg.genome_id]["style"] == 'dashed'):
            edge_flag = False
        elif g.has_edge(left_kmer_2, right_kmer, key = seg.genome_id) and g[left_kmer_2][right_kmer][seg.genome_id]["style"] == 'dashed':
            edge_flag = False
        elif g.has_edge(left_kmer, right_kmer_2, key = seg.genome_id) and g[left_kmer][right_kmer_2][seg.genome_id]["style"] == 'dashed':
            edge_flag = False
        elif g.has_edge(left_kmer_2, right_kmer_2, key = seg.genome_id)and g[left_kmer_2][right_kmer_2][seg.genome_id]["style"] == 'dashed':
            edge_flag = False
        elif edge_flag:
            if g.has_edge(left_kmer, right_kmer, key = seg.genome_id):
                g[left_kmer][right_kmer][seg.genome_id]["support"] =g[left_kmer][right_kmer][seg.genome_id]["support"]+ int(seg.coverage[0])
                g[left_kmer][right_kmer][seg.genome_id]["penwidth"] = 2
            else:
                g.add_edge(left_kmer, right_kmer, key=seg.genome_id, support=int(seg.coverage[0]), style='solid', penwidth = 1)
            label="L:{0} C:{1}".format(seg.length_bp, g[left_kmer][right_kmer][seg.genome_id]["support"])
            g[left_kmer][right_kmer][seg.genome_id]["label"] = label
            
    for double_bp in double_breaks:
        update_adjacencies(double_bp)
    
    ref_style = "dashed" if reference_adjacencies else "invis"
    for seg in genomicsegments:
        if seg.length_bp < max_genomic_len:
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
    print(key_to_color)

    return g, key_to_color


def add_legend(key_to_color, fout):
    fout.write("subgraph cluster_01 {\n\tlabel = \"Legend\";\n\tnode [shape=point]\n{\n\trank=same\n")
    for i, (key, color) in enumerate(key_to_color.items()):
        node_1 = "legend_{0}_1".format(i)
        node_2 = "legend_{0}_2".format(i)
        fout.write("\t{0} -- {1} [color={2}, label=\"{3}\"];\n".format(node_1, node_2, color, key))
    fout.write("}\n};\n")


def output_connected_components(graph, out_stream, target_genomes, control_genomes):
    components_list = []
    t= 0
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


def build_breakpoint_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies,
                           out_file, target_genomes, control_genomes):
    graph, key_to_color = build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies)

    with open(out_file, "w") as out_stream:
        out_stream.write("graph {\n")
        out_stream.write("edge [penwidth=2];\n")
        add_legend(key_to_color, out_stream)
        output_connected_components(graph, out_stream, target_genomes, control_genomes)
        out_stream.write("}\n")

    #tmp_graph = out_file + "_tmp"
    #nx.drawing.nx_pydot.write_dot(graph, tmp_graph)
    #os.remove(tmp_graph)

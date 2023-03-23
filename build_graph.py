#!/usr/bin/env python3

import argparse
import sys
import os
import networkx as nx
from collections import defaultdict
from copy import copy


#def neg_segment(segment):
#    if segment.startswith("+"):
#        return "-" + segment[1:]
#    else:
#        return "+" + segment[1:]

#COLORS = ["yellowgreen", "thistle", "peachpuff", "yellow", "khaki", "steelblue", "hotpink", "preu"]
COLORS = ["red", "green", "blue", "orange", "purple", "cyan", "hotpink", "preu"]
SEQUENCE_KEY = "__genomic"


def build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies):
    #kmer_size = 1
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]

    ### Adding adjacency edges
    for double_bp in double_breaks:
        left_kmer = node_to_id(double_bp.bp_1.unique_name())
        right_kmer = node_to_id(double_bp.bp_2.unique_name())

        if not g.has_node(left_kmer):
            g.add_node(left_kmer,label=double_bp.bp_1.fancy_name())
        if not g.has_node(right_kmer):
            g.add_node(right_kmer, label=double_bp.bp_2.fancy_name())

        if g.has_edge(left_kmer, right_kmer, key=double_bp.genome_id):
            g[left_kmer][right_kmer][double_bp.genome_id]["support"] += double_bp.supp
        else:
            #edge_weight = 2 if double_bp.genotype == "hom" else 1
            g.add_edge(left_kmer, right_kmer, key=double_bp.genome_id, support=double_bp.supp,
                       style=double_bp.edgestyle, adjacency=True, genotype=double_bp.genotype)

        read_support = g[left_kmer][right_kmer][double_bp.genome_id]["support"]
        g[left_kmer][right_kmer][double_bp.genome_id]["label"] = f"R:{read_support}"

    ### Linking complementary nodes + haplotype switches style
    compl_link_style = "dashed" if reference_adjacencies else "invis"
    for seg in genomicsegments:
        for coord in [seg.pos1, seg.pos2]:
            pos_bp, neg_bp = f"+{seg.ref_id}:{coord}", f"-{seg.ref_id}:{coord}"
            pos_node, neg_node = node_to_id(pos_bp), node_to_id(neg_bp)

            if not g.has_node(pos_node):
                g.add_node(pos_node, label=pos_bp)
            if not g.has_node(neg_node):
                g.add_node(neg_node, label=neg_bp)
            if coord not in hb_points[seg.ref_id]:
                g.add_edge(pos_node, neg_node, key=SEQUENCE_KEY, style=compl_link_style, adjacency=False, genotype="het")

    #helper function to add a new or update an existing genomic edge
    def _update_nx(left_node, left_label, left_break, right_node, right_label, right_break, genome_id, coverage):
        if g.has_edge(left_node, right_node, key=genome_id):
            g[left_node][right_node][genome_id]["support"] += coverage
            g[left_node][right_node][genome_id]["genotype"] = "hom"
        else:
            g.add_edge(left_node, right_node, key=genome_id,
                       support=coverage, style='solid', genotype="het", adjacency=False)

        support = g[left_node][right_node][seg.genome_id]["support"]
        g[left_node][right_node][seg.genome_id]["label"] = f"L:{seg.length_bp} C:{support}"

        g.nodes[left_node]["label"] = left_label
        if left_break:
            g.nodes[left_node]["style"] = "rounded, filled"
            g.nodes[left_node]["fillcolor"] = "grey"

        g.nodes[right_node]["label"] = right_label
        if right_break:
            g.nodes[right_node]["style"] = "rounded, filled"
            g.nodes[right_node]["fillcolor"] = "grey"

    ### Add genomic edges
    for seg_num, seg in enumerate(genomicsegments):
        if seg.haplotype == 0 and seg.coverage == 0:
            continue

        left_label, right_label = f"{seg.dir1}{seg.ref_id}:{seg.pos1}", f"{seg.dir2}{seg.ref_id}:{seg.pos2}"
        left_node, right_node = node_to_id(left_label), node_to_id(right_label)

        if seg.length_bp < max_genomic_len:
            left_break = seg.pos1 in hb_points[seg.ref_id]
            right_break = seg.pos2 in hb_points[seg.ref_id]
            if left_break:
                left_node = node_to_id(f"split_{seg_num}")
            if right_break:
                right_node = node_to_id(f"split_{seg_num}")
            _update_nx(left_node, left_label, left_break, right_node, right_label, right_break, seg.genome_id, seg.coverage)

        else:
            split_1, split_2 = node_to_id(f"split_{seg_num}_1"), node_to_id(f"split_{seg_num}_2")
            _update_nx(left_node, left_label, False, split_1, right_label, True, seg.genome_id, seg.coverage)
            _update_nx(split_2, left_label, True, right_node, right_label, False, seg.genome_id, seg.coverage)

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
            node_1 = f"legend_{i}_1"
            node_2 = f"legend_{i}_2"
            fout.write(f"\t{node_1} -> {node_2} [color={color}, label=\"{key}\", dir=none];\n")
        fout.write("}\n};\n")

    def _draw_components(components_list, out_stream):
        for subgr_num, (rank, cc, target_adj, control_adj, num_somatic) in enumerate(components_list):
            if num_somatic == 0:
                continue

            out_stream.write("subgraph cluster{0} {{\n".format(subgr_num))
            #Nodes properties
            for n in cc:
                properties = []
                for key, val in graph.nodes[n].items():
                    properties.append(f"{key}=\"{val}\"")
                prop_str = ",".join(properties)
                out_stream.write(f"{n} [{prop_str}];\n")

            #Edges
            for u, v, key in graph.edges(cc, keys=True):
                #switching nodes order to maintain genomic direction
                u_label = graph.nodes[u]["label"]
                v_label = graph.nodes[v]["label"]
                u_sign, (u_chr, u_pos) = u_label[0], u_label[1:].split(":")
                v_sign, (v_chr, v_pos) = v_label[0], v_label[1:].split(":")

                if ("INS" not in u_label) and ("INS" not in v_label):
                    if (v_chr, int(v_pos), v_sign) > (u_chr, int(u_pos), u_sign):
                        u, v = v, u
                elif "INS" in u_label and v_sign == "-":
                    u, v = v, u
                elif "INS" in v_label and u_sign == "+":
                    u, v = v, u

                #double edges for homozygous variants
                if graph[u][v][key]["genotype"] == "hom":
                    graph[u][v][key]["color"] = graph[u][v][key]["color"] + ":" + graph[u][v][key]["color"]

                #no arrows for adjacency edges + no weight in ranking
                properties = []
                if graph[u][v][key]["adjacency"] == True:
                    properties.append("dir=none")
                    properties.append("weight=0")

                for prop, val in graph[u][v][key].items():
                    properties.append(f"{prop}=\"{val}\"")
                prop_str = ",".join(properties)
                out_stream.write(f"{u} -> {v} [{prop_str}];\n")

            out_stream.write("}\n\n")

    with open(out_file, "w") as out_stream:
        out_stream.write("digraph {\n")
        #out_stream.write("edge [penwidth=2];\n")
        out_stream.write("node [shape=\"box\", style=\"rounded\"];\n")
        _add_legend(key_to_color, out_stream)
        _draw_components(connected_components, out_stream)
        out_stream.write("}\n")


def output_clusters_csv(graph, connected_components, out_file):
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\tadj_1\tadj_2\tgenome_ids\tread_support\tgenotype\n")
        for subgr_num, (rank, cc, target_adj, control_adj, _num_somatic) in enumerate(connected_components):
            all_adj = defaultdict(list)
            for u, v, key in graph.edges(cc, keys=True):
                if key != SEQUENCE_KEY and graph[u][v][key]["adjacency"] == True:
                    all_adj[(u, v)].append(key)

            for (u, v), keys in all_adj.items():
                keys.sort()
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
        num_somatic_adj = 0
        for u, v, key, data in graph.edges(cc, keys=True, data=True):
            if key in target_genomes and data["adjacency"]:
                target_adj.add((u, v))
            if key in control_genomes and data["adjacency"]:
                control_adj.add((u, v))

        for adj in target_adj:
            if adj not in control_adj:
                num_somatic_adj += 1

        rank = None
        if len(target_adj) > 0 and len(control_adj) == 0:
            rank = len(target_adj) + 100
        elif len(target_adj) == 0 and len(control_adj) == 0:
            rank = -100
        elif len(target_adj) == 0 and len(control_adj) > 0:
            rank = -len(control_adj)
        else:
            rank = len(target_adj) / len(control_adj)

        if len(target_adj) + len(control_adj):
            components_list.append((rank, cc, len(target_adj), len(control_adj), num_somatic_adj))

    components_list.sort(key=lambda p: p[0], reverse=True)

    return components_list


def build_breakpoint_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies,
                           target_genomes, control_genomes):
    graph, key_to_color = build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies)
    adj_clusters = cluster_adjacencies(graph, target_genomes, control_genomes)
    return graph, adj_clusters, key_to_color

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
            g.add_node(left_kmer, _coordinate=double_bp.bp_1.coord_tuple(), _terminal=False, _insertion=double_bp.bp_1.insertion_size)
        if not g.has_node(right_kmer):
            g.add_node(right_kmer, _coordinate=double_bp.bp_2.coord_tuple(), _terminal=False, _insertion=double_bp.bp_2.insertion_size)

        if not g.has_edge(left_kmer, right_kmer, key=double_bp.genome_id):
            g.add_edge(left_kmer, right_kmer, key=double_bp.genome_id, _support=double_bp.supp,
                       _type="adjacency", _genotype=double_bp.genotype)
        else:
            #we can have edges from the same genome, but different haplotypes
            g[left_kmer][right_kmer][double_bp.genome_id]["_support"] += double_bp.supp

    #helper function to add a new or update an existing genomic edge
    def _update_genomic(left_node, left_coords, left_terminal, right_node, right_coords, right_terminal, genome_id, coverage):
        if not g.has_node(left_node):
            g.add_node(left_node, _coordinate=left_coords, _terminal=left_terminal, _insertion=None)
        if not g.has_node(right_node):
            g.add_node(right_node, _coordinate=right_coords, _terminal=right_terminal, _insertion=None)

        if not g.has_edge(left_node, right_node, key=genome_id):
            g.add_edge(left_node, right_node, key=genome_id,
                       _support=coverage, _genotype="het", _type="genomic")
        else:
            #different haplotype
            g[left_node][right_node][genome_id]["_support"] += coverage
            g[left_node][right_node][genome_id]["_genotype"] = "hom"

    ### Add genomic edges
    for seg_num, seg in enumerate(genomicsegments):
        if seg.haplotype == 0 and seg.coverage == 0:
            continue

        left_label, right_label = seg.left_coord_str(), seg.right_coord_str()
        left_node, right_node = node_to_id(left_label), node_to_id(right_label)
        left_coord, right_coord = seg.left_coord_tuple(), seg.right_coord_tuple()

        if seg.length_bp < max_genomic_len:
            left_terminal = seg.pos1 in hb_points[seg.ref_id]
            right_terminal = seg.pos2 in hb_points[seg.ref_id]
            _update_genomic(left_node, left_coord, left_terminal, right_node, right_coord, right_terminal, seg.genome_id, seg.coverage)

        else:
            split_1, split_2 = node_to_id(right_label + "_split"), node_to_id(left_label + "_split")
            _update_genomic(left_node, left_coord, False, split_1, right_coord, True, seg.genome_id, seg.coverage)
            _update_genomic(split_2, left_coord, True, right_node, right_coord, False, seg.genome_id, seg.coverage)

    ### Linking complementary nodes + haplotype switches style
    #compl_link_style = "dashed" if reference_adjacencies else "invis"
    for seg in genomicsegments:
        for coord in [seg.pos1, seg.pos2]:
            pos_bp, neg_bp = f"+{seg.ref_id}:{coord}", f"-{seg.ref_id}:{coord}"
            pos_node, neg_node = node_to_id(pos_bp), node_to_id(neg_bp)

            if coord not in hb_points[seg.ref_id] and g.has_node(pos_node) and g.has_node(neg_node):
                g.add_edge(pos_node, neg_node, key=SEQUENCE_KEY, _type="complementary", _genotype="het")

    return g


def _node_to_str(node_dict, commas):
    ref, pos, sign = node_dict["_coordinate"]
    insertion = node_dict["_insertion"]
    if insertion is None:
        if not commas:
            return f"{sign}{ref}:{pos}"
        else:
            return f"{sign}{ref}:{pos:,}"
    else:
        return f"INS:{insertion}"


def _coord_cmp(node_data_1, node_data_2):
    if (node_data_1["_insertion"] is None) and (node_data_2["_insertion"] is None):
        return node_data_1["_coordinate"] < node_data_2["_coordinate"]
    elif node_data_1["_insertion"]:
        return node_data_2["_coordinate"][2] == "-"
    elif node_data_2["_insertion"]:
        return node_data_1["_coordinate"][2] == "+"


def output_clusters_graphvis(graph, connected_components, out_file):
    #node visual attributes
    for n in graph.nodes:
        graph.nodes[n]["label"] = _node_to_str(graph.nodes[n], commas=True)
        if graph.nodes[n]["_terminal"]:
            graph.nodes[n]["style"] = "rounded, filled"
            graph.nodes[n]["fillcolor"] = "grey"
        else:
            graph.nodes[n]["style"] = "rounded"

    #edges visual attributes
    for u, v, key in graph.edges:
        if graph[u][v][key]["_type"] == "adjacency":
            support = graph[u][v][key]["_support"]
            graph[u][v][key]["label"] = f"R:{support}"
            graph[u][v][key]["style"] = "dashed"
            graph[u][v][key]["dir"] = "none"

        elif graph[u][v][key]["_type"] == "genomic":
            coverage = graph[u][v][key]["_support"]
            #length = abs(graph.nodes[u]["_coordinate"][1] - graph.nodes[v]["_coordinate"][1])
            #graph[u][v][key]["label"] = f"L:{length}\\nC:{coverage}"
            graph[u][v][key]["label"] = f"C:{coverage}"

        elif graph[u][v][key]["_type"] == "complementary":
            graph[u][v][key]["style"] = "invis"

    #double edges for homozygous variants
    if graph[u][v][key]["_genotype"] == "hom":
        graph[u][v][key]["color"] = graph[u][v][key]["color"] + ":" + graph[u][v][key]["color"]

    #adding colors
    key_to_color = {}
    next_color = 0
    for u, v, _key in graph.edges:
        if _key != SEQUENCE_KEY:
            if _key not in key_to_color:
                key_to_color[_key] = COLORS[next_color]
                next_color = (next_color + 1) % len(COLORS)

            if graph[u][v][_key]["_genotype"] == "het":
                graph[u][v][_key]["color"] = key_to_color[_key]
            else:
                graph[u][v][_key]["color"] = key_to_color[_key] +":" + key_to_color[_key]

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
                if not _coord_cmp(graph.nodes[u], graph.nodes[v]):
                    u, v = v, u

                properties = []
                for prop, val in graph[u][v][key].items():
                    properties.append(f"{prop}=\"{val}\"")
                prop_str = ",".join(properties)
                out_stream.write(f"{u} -> {v} [{prop_str}];\n")

            out_stream.write("}\n\n")

    with open(out_file, "w") as out_stream:
        out_stream.write("digraph {\n")
        #out_stream.write("edge [penwidth=2];\n")
        out_stream.write("node [shape=\"box\"];\n")
        _add_legend(key_to_color, out_stream)
        _draw_components(connected_components, out_stream)
        out_stream.write("}\n")


def output_clusters_csv(graph, connected_components, out_file):
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\tadj_1\tadj_2\tgenome_ids\tread_support\tgenotype\n")
        for subgr_num, (rank, cc, target_adj, control_adj, _num_somatic) in enumerate(connected_components):
            all_adj = defaultdict(list)
            for u, v, key in graph.edges(cc, keys=True):
                if key != SEQUENCE_KEY and graph[u][v][key]["_type"] == "adjacency":
                    all_adj[(u, v)].append(key)

            for (u, v), keys in all_adj.items():
                if not _coord_cmp(graph.nodes[u], graph.nodes[v]):
                    u, v = v, u

                keys.sort()
                label_1 = _node_to_str(graph.nodes[u], commas=False)
                label_2 = _node_to_str(graph.nodes[v], commas=False)

                keys_text = ",".join(keys)
                read_support = ",".join([str(graph[u][v][k]["_support"]) for k in keys])
                genotypes = ",".join([graph[u][v][k]["_genotype"] for k in keys])
                fout.write(f"bga_{subgr_num}\t{label_1}\t{label_2}\t{keys_text}\t{read_support}\t{genotypes}\n")


def cluster_adjacencies(graph, target_genomes, control_genomes):
    components_list = []
    for cc in nx.connected_components(graph):
        target_adj = set()
        control_adj = set()
        num_somatic_adj = 0
        for u, v, key, data in graph.edges(cc, keys=True, data=True):
            if key in target_genomes and data["_type"] == "adjacency":
                target_adj.add((u, v))
            if key in control_genomes and data["_type"] == "adjacency":
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
    graph = build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies)
    adj_clusters = cluster_adjacencies(graph, target_genomes, control_genomes)
    return graph, adj_clusters

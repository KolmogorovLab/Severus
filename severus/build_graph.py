#!/usr/bin/env python3

import os
import networkx as nx
from collections import defaultdict
from itertools import combinations
import logging

from severus.breakpoint_finder import get_genomic_segments, cluster_db, add_inbetween_ins, output_breaks
from severus.vcf_output import write_to_vcf

logger = logging.getLogger()


#COLORS = ["yellowgreen", "thistle", "peachpuff", "yellow", "khaki", "steelblue", "hotpink", "preu"]
COLORS = ["#189BA0", "#830042", "#B2C971", "#8470FF", "#1B80B3", "#FF7A33", "#B35900", "#006400"]
SEQUENCE_KEY = "__genomic"

def filter_small_sv(double_breaks):
    db_to_remove = []
    MIN_SV_SIZE = 5000
    for db in double_breaks:
        if db.bp_1.ref_id == db.bp_2.ref_id and not db.sv_type == 'foldback' and abs(db.bp_1.position - db.bp_2.position) < MIN_SV_SIZE:
            db_to_remove.append(db)
        elif db.bp_2.is_insertion and db.bp_2.insertion_size < MIN_SV_SIZE:
            db_to_remove.append(db)
            
    if db_to_remove:
        for db in list(set(db_to_remove)):
            double_breaks.remove(db) 
                    
def build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies, filter_small_svs):
    if filter_small_svs:
        filter_small_sv(double_breaks)
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    
    def node_to_id(node_str):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        return node_ids[node_str]
    
    #helper function to add a new or update an existing genomic edge
    def _update_genomic(left_node, left_coord, left_terminal, right_node, right_coord, right_terminal, genome_id, coverage, haplotype):
        if not g.has_node(left_node):
            g.add_node(left_node, _coordinate=left_coord[:2], _loose_end = False, _terminal=left_terminal, _insertion=None, _contig=None, _dir = left_coord[2], _cluster_id = None)
        if not g.has_node(right_node):
            g.add_node(right_node, _coordinate=right_coord[:2], _loose_end = False, _terminal=right_terminal, _insertion=None,_contig=None, _dir = right_coord[2], _cluster_id = None)
        
        if g.has_edge(left_node, right_node, key=genome_id) and g[left_node][right_node][genome_id]['_type'] == 'genomic':
            g[left_node][right_node][genome_id]["_support"] += coverage
            g[left_node][right_node][genome_id]["_haplotype"] += str(haplotype)
            if not haplotype == 0:
                g[left_node][right_node][genome_id]["_genotype"] = "hom"
       
        elif not g.has_edge(left_node, right_node, key=genome_id):
            g.add_edge(left_node, right_node, key=genome_id,
                       _support=coverage, _genotype="het", _type="genomic", _haplotype = str(haplotype))
            
    ### Adding adjacency edges
    for double_bp in double_breaks:
        left_kmer = node_to_id(double_bp.bp_1.unique_name())
        right_kmer = node_to_id(double_bp.bp_2.unique_name())
        
        if not g.has_node(left_kmer):
            g.add_node(left_kmer, _coordinate=double_bp.bp_1.coord_tuple()[:2], _loose_end = False, _terminal=False, _insertion=double_bp.bp_1.insertion_size,_contig = double_bp.bp_1.contig_id,_dir = double_bp.bp_1.dir_1, _cluster_id = double_bp.subgraph_id)
        if not g.has_node(right_kmer):
            g.add_node(right_kmer, _coordinate=double_bp.bp_2.coord_tuple()[:2], _loose_end = double_bp.bp_2.loose_end_id, _terminal=False, _insertion=double_bp.bp_2.insertion_size, _contig = double_bp.bp_2.contig_id,_dir = double_bp.bp_2.dir_1, _cluster_id = double_bp.subgraph_id)
        
        if not g.has_edge(left_kmer, right_kmer, key=double_bp.genome_id):
            g.add_edge(left_kmer, right_kmer, key=double_bp.genome_id, _support=double_bp.supp,
                       _type="adjacency", _genotype=double_bp.genotype)
        else:
            #we can have edges from the same genome, but different haplotypes                                                 
            g[left_kmer][right_kmer][double_bp.genome_id]["_support"] += double_bp.supp

    ### Add genomic edges
    for seg_num, seg in enumerate(genomicsegments):                                                                                                                                                                                                                                                                                          
        if seg.haplotype == 0 and seg.coverage == 0:
            continue
        
        left_label, right_label = seg.left_coord_str(), seg.right_coord_str()
        left_node, right_node = node_to_id(left_label), node_to_id(right_label)
        left_coord, right_coord = seg.left_coord_tuple(), seg.right_coord_tuple()
        
        same_clust = False
        if g.has_node(left_node) and g.has_node(right_node):
            same_clust = True if g.nodes[left_node]['_cluster_id'] and g.nodes[right_node]['_cluster_id'] and g.nodes[left_node]['_cluster_id'] == g.nodes[right_node]['_cluster_id'] else False
            
        if seg.length_bp < max_genomic_len or same_clust:
            left_terminal = seg.pos1 in hb_points[seg.ref_id]
            if left_terminal:
                left_node = node_to_id(f"hb_left_{left_label}")
            right_terminal = seg.pos2 in hb_points[seg.ref_id]
            if right_terminal:
                right_node = node_to_id(f"hb_right_{right_label}")
            _update_genomic(left_node, left_coord, left_terminal, right_node, right_coord, right_terminal, seg.genome_id, seg.coverage, seg.haplotype)
            
        else:
            split_1, split_2 = node_to_id(f"split_1_{left_label}"), node_to_id(f"split_2_{right_label}")
            _update_genomic(left_node, left_coord, False, split_1, right_coord, True, seg.genome_id, seg.coverage, seg.haplotype)
            _update_genomic(split_2, left_coord, True, right_node, right_coord, False, seg.genome_id, seg.coverage, seg.haplotype)
        
        
    return g


def _node_to_str(node_dict, commas):
    ref, pos= node_dict["_coordinate"]
    insertion = node_dict["_insertion"]
    if insertion is None:
        if not commas:
            return f"{ref}:{pos}"
        else:
            return f"{ref}:{pos:,}"
    else:
        return f"INS:{insertion}"

def _node_to_str_output(node_dict, commas):
    ref, pos= node_dict["_coordinate"]
    sign = '-' if node_dict["_dir"] == -1 else '+'
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
    if node_data_1["_contig"] or node_data_2["_contig"]:
        return node_data_1["_coordinate"] < node_data_2["_coordinate"]
    if node_data_1["_insertion"] and not node_data_1["_contig"]:
        return True

    
def output_clusters_graphvis(graph, connected_components, out_file):
    
    def _add_legend(key_to_color, fout):
        fout.write("subgraph cluster_01 {\n\tlabel = \"Legend\";\n\tnode [shape=point]\n{\n")
        for i, (key, color) in enumerate(key_to_color.items()):
            node_1 = f"legend_{i}_1"
            node_2 = f"legend_{i}_2"
            fout.write(f"\t{node_1} -> {node_2} [color=\"{color}\", label=\"{key}\", dir=none];\n")
        fout.write("}\n};\n")
    
    def _add_dir(dir_1):
        return 'w' if dir_1 == -1 else 'e'

    def _draw_components(graph, components_list, out_stream):
        subgr_num = len(components_list)
        for (rank, cc, target_adj, control_adj, num_somatic) in reversed(components_list):
            subgr_num -= 1
            if num_somatic == 0:
                continue
            
            out_stream.write("subgraph cluster{0} {{\n".format(subgr_num))
            out_stream.write("label = \"subgraph cluster {0}\"\nlabeljust=l\n".format(subgr_num))
            subgr_chr = defaultdict(list)
            
            #Nodes properties
            for n in cc:
                properties = []
                t = 0
                for key, val in graph.nodes[n].items():
                    properties.append(f"{key}=\"{val}\"")
                prop_str = ",".join(properties)
                if not graph.nodes[n]['_insertion']:
                    subgr_chr[graph.nodes[n]['_coordinate'][0]].append((n,graph.nodes[n]['_coordinate'][1]))
                out_stream.write(f"{n} [{prop_str}];\n")
                
            for val in subgr_chr.values():
                if len(val) == 1:
                    continue
                tt = ' -> '.join([str(x[0]) for x in sorted(val, key = lambda x:x[1])])
                out_stream.write(f"{tt} [style = invis ];\n")
                
            #Edges    
            for u, v, key in graph.edges(cc, keys=True):
                #switching nodes order to maintain genomic direction
                if not _coord_cmp(graph.nodes[u], graph.nodes[v]):
                    u, v = v, u
                    
                properties = []
                for prop, val in graph[u][v][key].items():
                    properties.append(f"{prop}=\"{val}\"")
                prop_str = ",".join(properties)
                if graph[u][v][key]['_type'] == 'genomic' or graph.nodes[u]['_contig'] or graph.nodes[v]['_contig']:
                    dir_1, dir_2 = 'e', 'w'
                elif graph.nodes[u]['_insertion'] and not graph.nodes[u]['_contig'] and not graph.nodes[v]['_contig']:
                    dir_1, dir_2 = 'n', 's'
                    out_stream.write("subgraph cluster{0}_{1} {{\ncolor=invis;\nlabel=\"\";\n{{rank=same; {2}; {3} }};\n}}".format(subgr_num, t, u, v))
                    t+=1
                else:
                    dir_1, dir_2 = _add_dir(graph.nodes[u]['_dir']),_add_dir(graph.nodes[v]['_dir'])
                out_stream.write(f"{u}:{dir_1} -> {v}:{dir_2} [{prop_str}];\n")
                
            out_stream.write("}\n\n")  
    
    #node visual attributes
    for n in graph.nodes:
        graph.nodes[n]["label"] = _node_to_str(graph.nodes[n], commas=True)
        if graph.nodes[n]["_terminal"]:
            graph.nodes[n]["style"] = "filled"
            graph.nodes[n]["fillcolor"] = "grey"
        if graph.nodes[n]["_loose_end"]:
            graph.nodes[n]["shape"] = 'point'

    #edges visual attributes
    k_THR = 2000
    for u, v, key in graph.edges:
        if graph[u][v][key]["_type"] == "adjacency":
            support = graph[u][v][key]["_support"]
            graph[u][v][key]["label"] = f"R:{support}"
            graph[u][v][key]["style"] = "dashed"
            graph[u][v][key]["dir"] = "both"
            graph[u][v][key]["arrowhead"] = 'box' if graph.nodes[u]['_dir'] == -1 else 'normal'
            graph[u][v][key]["arrowtail"] = 'box' if graph.nodes[v]['_dir'] == -1 else 'normal'
        elif graph[u][v][key]["_type"] == "genomic":
            coverage = graph[u][v][key]["_support"]
            len_g = abs(graph.nodes[u]["_coordinate"][1] - graph.nodes[v]["_coordinate"][1])
            if len_g > k_THR:
                length_1k = len_g / 1000
                graph[u][v][key]["label"] = f"C:{coverage}\\nL:{length_1k:.1f}kb"
            else:
                graph[u][v][key]["label"] = f"C:{coverage}\\nL:{len_g}bp"
        elif graph[u][v][key]["_type"] == "complementary":
            graph[u][v][key]["style"] = "invis"

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
                graph[u][v][_key]["color"] = key_to_color[_key] + ":" + key_to_color[_key]

    with open(out_file, "w") as out_stream:
        out_stream.write("digraph {\n")
        out_stream.write("node [shape=\"box\"];\n rankdir = \"LR\"")
        _draw_components(graph, connected_components, out_stream)
        _add_legend(key_to_color, out_stream)
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

                #only output insertions once in text
                if graph.nodes[v]["_insertion"]:
                    continue

                keys.sort()
                label_1 = _node_to_str_output(graph.nodes[u], commas=False)
                label_2 = _node_to_str_output(graph.nodes[v], commas=False)

                keys_text = ",".join(keys)
                read_support = ",".join([str(graph[u][v][k]["_support"]) for k in keys])
                genotypes = ",".join([graph[u][v][k]["_genotype"] for k in keys])
                fout.write(f"severus_{subgr_num}\t{label_1}\t{label_2}\t{keys_text}\t{read_support}\t{genotypes}\n")

def cc_to_label(graph, connected_components):
    id_to_cc = defaultdict(int)
    for subgr_num, (rank, cc, _, _, _) in enumerate(connected_components):
        all_adj = defaultdict(list)
        for u, v, key in graph.edges(cc, keys=True):
            if key != SEQUENCE_KEY and graph[u][v][key]["_type"] == "adjacency":
                all_adj[(u, v)].append(key)
        for (u, v), keys in all_adj.items():
            if graph.nodes[u]["_insertion"]:
                continue
            keys.sort()
            label_1 = _node_to_str_output(graph.nodes[u], commas=False)
            label_2 = _node_to_str_output(graph.nodes[v], commas=False)
            id_to_cc[label_1 + '|' + label_2] = subgr_num
    return id_to_cc

def cluster_adjacencies(graph, target_genomes, control_genomes):
    components_list = []
    for cc in nx.connected_components(graph):
        target_adj = set()
        control_adj = set()
        num_somatic_adj = 0
        for u, v, key, data in graph.edges(cc, keys=True, data=True):
            if graph.nodes[u]['_loose_end'] or graph.nodes[v]['_loose_end']:
                continue
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
                           target_genomes, control_genomes, filter_small_svs):
    graph = build_graph(double_breaks, genomicsegments, hb_points, max_genomic_len, reference_adjacencies, filter_small_svs)
    adj_clusters = cluster_adjacencies(graph, target_genomes, control_genomes)
    return graph, adj_clusters
            
            
def output_graphs(db_list, coverage_histograms, thread_pool, target_genomes, control_genomes, ref_lengths, args):
    keys = ['germline', 'somatic']
    for key in keys:
        double_breaks = db_list[key]
        if key == 'germline' and args.only_somatic:
            continue
        
        sub_fol = "all_SVs" if key == 'germline' else 'somatic_SVs'
        out_folder = os.path.join(args.out_dir, sub_fol)
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
        all_ids = target_genomes + control_genomes if key == 'germline' else target_genomes   
            
        logger.info(f"Preparing outputs for {sub_fol}")
        
        out_breakpoint_graph = os.path.join(out_folder, "breakpoint_graph.gv")
        out_clustered_breakpoints = os.path.join(out_folder, "breakpoint_clusters.csv")
        
        logger.info("\tComputing segment coverage")
        genomic_segments, hb_points = get_genomic_segments(double_breaks, coverage_histograms, thread_pool, args.phase_vcf)
        genomic_segments.sort(key = lambda seg:seg.haplotype, reverse = True )
        
        if args.inbetween_ins and key == 'germline':
            double_breaks2 = add_inbetween_ins(double_breaks)
        else:
            double_breaks2 = double_breaks
        
        logger.info("\tPreparing graph")
        graph, adj_clusters = build_breakpoint_graph(double_breaks2, genomic_segments, hb_points, args.max_genomic_len,
                                                                    args.reference_adjacencies, target_genomes, control_genomes, args.filter_small_svs)
        output_clusters_graphvis(graph, adj_clusters, out_breakpoint_graph)
        output_clusters_csv(graph, adj_clusters, out_clustered_breakpoints)
        
        logger.info("\tWriting vcf")
        id_to_cc = cc_to_label(graph, adj_clusters)
        write_to_vcf(double_breaks, all_ids, id_to_cc, out_folder, key, ref_lengths, args.only_somatic, args.no_ins)
        
    
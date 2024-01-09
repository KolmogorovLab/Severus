#!/usr/bin/env python3

import os
import networkx as nx
from collections import defaultdict
import logging
import numpy as np

from severus.breakpoint_finder import get_genomic_segments, cluster_inversions, add_inbetween_ins
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
                    
def build_graph(genomic_segments, phasing_segments, adj_segments, max_genomic_len, reference_adjacencies):
    g = nx.MultiGraph()

    node_ids = {}
    id_to_kmers = {}
    db_to_cl = defaultdict(list)
    
    def node_to_id(node_str,db):
        if not node_str in node_ids:
            node_ids[node_str] = len(node_ids) + 1
            id_to_kmers[node_ids[node_str]] = node_str
        if db:
            db_to_cl[node_ids[node_str]].append(db)
        return node_ids[node_str]
    
    for db, gs in genomic_segments.items():
        if len(gs) < 3:
            left_kmer = node_to_id(gs[0].full_name(),db)
            right_kmer = node_to_id(gs[-1].full_name(),db)
            if not g.has_node(left_kmer):
                g.add_node(left_kmer, _coordinate=gs[0].full_name(), _coordinate_tuple = (gs[0].ref_id, gs[0].pos1, gs[0].pos2), _loose_end = False, _terminal=False, _insertion=db.bp_1.insertion_size,
                           _phase_switch = False, _coverage = gs[0].coverage, _length = gs[0].length_bp)
            if not g.has_node(right_kmer):
                g.add_node(right_kmer, _coordinate=gs[-1].full_name(), _coordinate_tuple = (gs[-1].ref_id, gs[-1].pos1, gs[-1].pos2), _loose_end = False, _terminal=False, _insertion=db.bp_2.insertion_size,
                           _phase_switch = False, _coverage = gs[-1].coverage, _length = gs[-1].length_bp)
            if not g.has_edge(left_kmer, right_kmer, key=db.genome_id):
                st1 = 'w' if db.bp_1.position in [gs[0].pos1, gs[-1].pos1] else 'e'
                st2 = 'w' if db.bp_2.position in [gs[0].pos1, gs[-1].pos1] else 'e'
                g.add_edge(left_kmer, right_kmer, key=db.genome_id, _support=db.supp,
                           _type="adjacency", _genotype=db.genotype, _dir = (st1, st2))
            else:                                      
                g[left_kmer][right_kmer][db.genome_id]["_support"] += db.supp
        if len(gs) == 3:
            ins_kmer = node_to_id(gs[0].full_name(),db)
            left_kmer = node_to_id(gs[1].full_name(),db)
            right_kmer = node_to_id(gs[2].full_name(),db)
            if not g.has_node(left_kmer):
                g.add_node(left_kmer, _coordinate=gs[1].full_name(), _coordinate_tuple = (gs[1].ref_id, gs[1].pos1, gs[1].pos2-1), _loose_end = False, _terminal=False, _insertion='',
                           _phase_switch = False, _coverage = gs[1].coverage, _length = gs[1].length_bp)
            if not g.has_node(right_kmer):
                g.add_node(right_kmer, _coordinate=gs[2].full_name(), _coordinate_tuple = (gs[2].ref_id, gs[2].pos1+1, gs[2].pos2), _loose_end = False, _terminal=False, _insertion='',
                           _phase_switch = False, _coverage = gs[2].coverage, _length = gs[2].length_bp)
            if not g.has_node(ins_kmer):
                g.add_node(ins_kmer, _coordinate=gs[0].ins_label(), _coordinate_tuple = (gs[1].ref_id, gs[1].pos2, gs[1].pos2), _loose_end = False, _terminal=False, _insertion='',
                           _phase_switch = False, _coverage = None, _length = None)
            if not g.has_edge(left_kmer, ins_kmer, key=db.genome_id):
                g.add_edge(left_kmer, ins_kmer, key=db.genome_id, _support=db.supp,
                           _type="adjacency", _genotype=db.genotype, _dir = ('e','w'))
            if not g.has_edge(ins_kmer, right_kmer, key=db.genome_id):
                g.add_edge(ins_kmer, right_kmer, key=db.genome_id, _support=db.supp,
                           _type="adjacency", _genotype=db.genotype, _dir = ('e','w'))
    
    for sw_p, gs_ls in phasing_segments.items():
        sw_p_kmer = node_to_id(sw_p[0] + ':' + str(sw_p[1]), '')
        if not g.has_node(sw_p_kmer):
            g.add_node(sw_p_kmer, _coordinate=sw_p, _coordinate_tuple = (sw_p[0],sw_p[1], sw_p[1]), _loose_end = False, _terminal=False, _insertion = None,
                       _phase_switch = True, _coverage = None, _length = None)
        for gs in gs_ls:
            left_kmer = node_to_id(gs.full_name(), '')
            if not g.has_edge(left_kmer, sw_p_kmer, key=gs.genome_id):
                st1 = ('w', 'e') if gs.pos1 == sw_p[1] else ('e', 'w')
                g.add_edge(left_kmer, sw_p_kmer, key=gs.genome_id, _support=0,
                           _type="block_sw", _genotype='', _dir = st1)
                
    for (gs1,gs2) in adj_segments:
        left_kmer = node_to_id(gs1.full_name(), '')
        right_kmer = node_to_id(gs2.full_name(), '')
        if g.has_node(left_kmer) and g.has_node(right_kmer):
            if not g.has_edge(left_kmer, right_kmer, key=gs1.genome_id):
                st1 = ('e', 'w')
                g.add_edge(left_kmer, right_kmer, key=gs1.genome_id, _support=0,
                           _type="ref_adj", _genotype= '', _dir = st1)
                
    return g, db_to_cl


def _node_to_str(node_dict):
    k_THR = 2000
    lab = node_dict["_coordinate"]
    if 'INS' in lab or node_dict["_phase_switch"]:
        return node_dict["_coordinate"]
    cov = node_dict["_coverage"]
    len_g = node_dict["_length"]
    if len_g > k_THR:
        length_1k = len_g / 1000
        label = f"C:{cov}, L:{length_1k:.1f}kb"
    else:
        label = f"C:{cov}, L:{len_g}bp" 
    return f"{lab}\\n{label}"

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
    return node_data_1['_coordinate_tuple'][1] < node_data_2['_coordinate_tuple'][1]
    
def output_clusters_graphvis(graph, connected_components, out_file):
    
    def _add_legend(key_to_color, fout):
        out_stream.write("digraph cluster_01 {\n")
        out_stream.write("node [shape=\"point\"];\n rankdir = \"LR\"")
        out_stream.write("subgraph cluster_01 {\n")
        out_stream.write("pencolor=transparent; \nlabel = \"Legend\"\nlabeljust=l\n")
        for i, (key, color) in enumerate(key_to_color.items()):
            node_1 = f"legend_{i}_1"
            node_2 = f"legend_{i}_2"
            fout.write(f"\t{node_1} -> {node_2} [color=\"{color}\", label=\"{key}\", dir=none];\n")
        fout.write("}\n}\n")

    def _draw_components(graph, components_list, out_stream):
        for subgr_num, (_,_,_, cc,_) in enumerate(components_list):
            
            out_stream.write("digraph cluster{0} {{\n".format(subgr_num))
            out_stream.write("node [shape=\"box\"];\n rankdir = \"LR\"")
            out_stream.write("subgraph cluster{0} {{\n".format(subgr_num))
            out_stream.write("pencolor=transparent; \nlabel = \"subgraph cluster {0}\"\nlabeljust=l\n".format(subgr_num))
            subgr_chr = defaultdict(list)
            
            #Nodes properties
            for n in cc:
                properties = []
                for key, val in graph.nodes[n].items():
                    properties.append(f"{key}=\"{val}\"")
                prop_str = ",".join(properties)
                if not graph.nodes[n]['_insertion']:
                    subgr_chr[graph.nodes[n]['_coordinate_tuple'][0]].append((n,graph.nodes[n]['_coordinate_tuple'][1]))
                out_stream.write(f"{n} [{prop_str}];\n")
                
            for val in subgr_chr.values():
                if len(val) == 1:
                    continue
                tt = ' -> '.join([str(x[0]) for x in sorted(val, key = lambda x:x[1])])
                out_stream.write(f"{tt} [style = invis ];\n")
                
            #Edges    
            for u, v, key, value in graph.edges(cc, keys=True, data = True):
                
                if value['_type'] == "cluster_conn":
                    continue
                
                if not _coord_cmp(graph.nodes[u], graph.nodes[v]):
                    u, v = v, u
                properties = []
                for prop, val in graph[u][v][key].items():
                    properties.append(f"{prop}=\"{val}\"")
                prop_str = ",".join(properties)
                dir_1, dir_2 = value['_dir']
                out_stream.write(f"{u}:{dir_1} -> {v}:{dir_2} [{prop_str}, dir=none];\n")
                
            out_stream.write("}\n}\n")  
    
    #node visual attributes
    for n in graph.nodes:
        graph.nodes[n]["label"] = _node_to_str(graph.nodes[n])
        if graph.nodes[n]["_terminal"]:
            graph.nodes[n]["style"] = "filled"
            graph.nodes[n]["fillcolor"] = "grey"
        if graph.nodes[n]['_phase_switch']:
            graph.nodes[n]["style"] = "filled"
            graph.nodes[n]["fillcolor"] = "grey"
            graph.nodes[n]["shape"] = 'point'

    #edges visual attributes
    for u, v, key in graph.edges:
        if graph[u][v][key]["_type"] == "adjacency":
            support = graph[u][v][key]["_support"]
            graph[u][v][key]["label"] = f"R:{support}"
            graph[u][v][key]["style"] = "dashed"
        elif graph[u][v][key]["_type"] == "ref_adj" or graph[u][v][key]["_type"] == "block_sw":
            graph[u][v][key]["style"] = "dashed"
            graph[u][v][key]["color"] = "grey"
            
    #adding colors
    key_to_color = {}
    next_color = 0
    for u, v, _key in graph.edges:
        if graph[u][v][_key]["_type"] == "adjacency":
            if _key not in key_to_color:
                key_to_color[_key] = COLORS[next_color]
                next_color = (next_color + 1) % len(COLORS)
                
            if graph[u][v][_key]["_genotype"] == "het":
                graph[u][v][_key]["color"] = key_to_color[_key]
            else:
                graph[u][v][_key]["color"] = key_to_color[_key] + ":" + key_to_color[_key]

    with open(out_file, "w") as out_stream:
        _add_legend(key_to_color, out_stream)
        _draw_components(graph, connected_components, out_stream)

def output_clusters_csv(db_to_cl, double_breaks, connected_components, out_file):
    
    for subgr_num, (_,_,_,cc,_) in enumerate(connected_components):
        
        for node_id in cc:
            for db in db_to_cl[node_id]:
                db.cluster_id = subgr_num
                
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\tadj_1\tadj_2\tgenome_ids\tread_support\tgenotype\n")            
        double_breaks.sort(key=lambda b:(b.cluster_id, b.genome_id))
        clusters = defaultdict(list) 
        
        for br in double_breaks:
            clusters[br.to_string()].append(br)
            
        for cl in clusters.values():
            supp_ls = defaultdict(int)
            gen_ls = defaultdict(str)
            
            for db in cl:
                supp_ls[db.genome_id] += db.supp
                gen_ls[db.genome_id] = db.genotype
                
            label_1 = db.to_string()
            keys_text = ",".join(list(gen_ls.keys()))
            read_support = ",".join([str(k) for k in supp_ls.values()])
            genotypes = ",".join(list(gen_ls.values()))
            fout.write(f"severus_{db.cluster_id}\t{label_1}\t{keys_text}\t{read_support}\t{genotypes}\n")

def output_clusters_info(connected_components, out_file):
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\type\tSVs\tSV_count\tgenome_ids\t\n") 
        for subgr_num, (sv_events,_type,_,_,genome) in enumerate(connected_components):
            nsv = sum([k for k in sv_events.values()])
            svs = ','.join([ f'{key}:{val} ' for key, val in sv_events.items() if val > 0])
            genome_ids = ','.join(genome)
            fout.write(f"severus_{subgr_num}\t{_type}\t{svs}\t{nsv}\t{genome_ids}\t\n")

def cluster_adjacencies(graph, db_to_cl, target_genomes, control_genomes):
    DIFF_THR = 10000
    INDEL_THR = 1.5
    components_list = []
    connected = defaultdict(list)
    
    for node, cl in db_to_cl.items():
        if cl[0].cluster_id:
            connected[cl[0].cluster_id].append(node)
    
    for cl in connected.values():
        for (a,b) in zip(cl[:-1],cl[1:]):
            db = db_to_cl[a][0]
            if not graph.has_edge(a,b, key=db.genome_id):
                graph.add_edge(a,b, key=db.genome_id, _support=0,
                           _type="cluster_conn", _genotype=db.genotype, _dir = ('e','w'))

    rank_ls = defaultdict(int)
    for i,cc in enumerate(nx.connected_components(graph)):
        sv_type = defaultdict(int)
        chr_list = defaultdict(list)
        sv_len = defaultdict(list)
        db_ls = []
        for node in cc:
            if graph.nodes[node]['_phase_switch']:
                continue
            db_ls.append(db_to_cl[node][0])
        db_ls = list(set(db_ls))
        genome = list(set([db.genome_id for db in db_ls]))
        if len(db_ls) == 1:
            db = db_ls[0]
            _cp = db.vcf_sv_type if not db.sv_type else db.sv_type
            sv_type[_cp] +=1
            _type = 'simple'
            if not db.bp_1.ref_id == db.bp_2.ref_id:
                _type = 'tra'
                rank = 1
            elif db.sv_type:
                rank = 2
            else:
                rank = 0
        else:
            _type = 'complex'
            rank = len(db_ls) + 2
            for db in db_ls:
                _cp = db.vcf_sv_type if not db.sv_type else db.sv_type
                sv_type[_cp] +=1
                chr_list[db.bp_1.ref_id].append(db.bp_1.position)
                if not db.bp_2.is_insertion:
                    chr_list[db.bp_2.ref_id].append(db.bp_2.position)
                if db.length > 0:
                    sv_len[db.bp_1.ref_id].append(db.length)
            if sv_type['DEL'] + sv_type['INS'] + sv_type['INV'] == len(db_ls):
                seq = list(chr_list.keys())[0]
                pos = chr_list[seq]
                pos.sort()
                diff = [a-b for (a,b) in zip(pos[1:],pos[:-1])]
                if np.median(diff) > np.median(sv_len[seq]) * INDEL_THR:
                    _type = 'simple'
                    rank = 3
        components_list.append((sv_type, _type, rank, cc, genome))
        rank_ls[rank] +=1
    
    components_list.sort(key=lambda p: p[2], reverse=True)
    
    return components_list


def build_breakpoint_graph(genomic_segments, phasing_segments, adj_segments, max_genomic_len, reference_adjacencies,
                           target_genomes, control_genomes):
    graph, db_to_cl = build_graph(genomic_segments, phasing_segments, adj_segments, max_genomic_len, reference_adjacencies)
    adj_clusters = cluster_adjacencies(graph, db_to_cl, target_genomes, control_genomes)
    return graph, adj_clusters, db_to_cl
            
            
def output_graphs(db_list, coverage_histograms, thread_pool, target_genomes, control_genomes, ref_lengths, args):
    keys = ['germline', 'somatic']
    for key in keys:
        if not key in list(db_list.keys()):
            continue
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
        out_cluster_list = os.path.join(out_folder, "breakpoint_clusters_list.tsv")
        
        logger.info("\tComputing segment coverage")
        if args.inbetween_ins and key == 'germline':
            double_breaks2 = add_inbetween_ins(double_breaks)
        else:
            double_breaks2 = double_breaks
            
        (genomic_segments, phasing_segments, adj_segments) = get_genomic_segments(double_breaks2, coverage_histograms, args.phase_vcf, key, ref_lengths, args.min_ref_flank, args.max_genomic_len)
        
        cluster_inversions(double_breaks2, coverage_histograms, args.min_sv_size)
        
        logger.info("\tPreparing graph")
        graph, adj_clusters, db_to_cl = build_breakpoint_graph(genomic_segments, phasing_segments, adj_segments, args.max_genomic_len,
                                                                    args.reference_adjacencies, target_genomes, control_genomes)
        output_clusters_graphvis(graph, adj_clusters, out_breakpoint_graph)
        output_clusters_csv(db_to_cl, double_breaks, adj_clusters, out_clustered_breakpoints)
        
        output_clusters_info(adj_clusters, out_cluster_list)
        
        logger.info("\tWriting vcf")
        write_to_vcf(double_breaks, all_ids, out_folder, key, ref_lengths, args.only_somatic, args.no_ins)
        
            
        
    
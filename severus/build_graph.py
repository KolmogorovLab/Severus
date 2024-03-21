#!/usr/bin/env python3

import os
import networkx as nx
from collections import defaultdict
import logging
import numpy as np
import plotly
import plotly.graph_objects as go

from severus.breakpoint_finder import get_genomic_segments, cluster_indels
from severus.vcf_output import write_to_vcf

logger = logging.getLogger()


#COLORS = ["yellowgreen", "thistle", "peachpuff", "yellow", "khaki", "steelblue", "hotpink", "preu"]
COLORS = ["#189BA0", "#830042", "#B2C971", "#8470FF", "#1B80B3", "#FF7A33", "#B35900", "#006400"]
SEQUENCE_KEY = "__genomic"

                    
def build_graph(genomic_segments,  adj_segments):

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
           
    ### graph build start
    for db, gs in genomic_segments.items():
        rn = [g for g in gs if db.bp_1.position in [g.pos1, g.pos2]]
        ln = [g for g in gs if db.bp_2.position in [g.pos1, g.pos2]]
        for r_node in rn:
            for l_node in ln:
                left_kmer = node_to_id(r_node.full_name(),db)
                right_kmer = node_to_id(l_node.full_name(),db)
                if not g.has_node(left_kmer):
                    g.add_node(left_kmer, _coordinate=r_node.full_name(), _coordinate_tuple = (r_node.ref_id, r_node.pos1, r_node.pos2), _loose_end = False, _terminal=False, _insertion=db.bp_1.insertion_size,
                               _phase_switch = False, _coverage = r_node.total_coverage, _hcoverage = r_node.coverage, _haplotype = [r_node.haplotype], _length = r_node.length_bp)
                else:
                    g.nodes[left_kmer]['_haplotype'].append(r_node.haplotype)
                    
                if not g.has_node(right_kmer):
                    g.add_node(right_kmer, _coordinate=l_node.full_name(), _coordinate_tuple = (l_node.ref_id, l_node.pos1, l_node.pos2), _loose_end = False, _terminal=False, _insertion=db.bp_2.insertion_size,
                               _phase_switch = False, _coverage = l_node.total_coverage,_hcoverage = l_node.coverage, _haplotype = [l_node.haplotype] , _length = l_node.length_bp)
                else:
                    g.nodes[right_kmer]['_haplotype'].append(l_node.haplotype)
                if not g.has_edge(left_kmer, right_kmer, key=db.genome_id):
                    st1 = 'w' if db.bp_1.position in [r_node.pos1, l_node.pos1] else 'e'
                    st2 = 'w' if db.bp_2.position in [r_node.pos1, l_node.pos1] else 'e'
                    g.add_edge(left_kmer, right_kmer, key=db.genome_id, _support=db.supp,
                               _type="adjacency", _genotype=db.genotype, _dir = (st1, st2))
                else:                                      
                    g[left_kmer][right_kmer][db.genome_id]["_support"] += db.supp
                
    for (gs1,gs2) in adj_segments:
        left_kmer = node_to_id(gs1.full_name(), '')
        right_kmer = node_to_id(gs2.full_name(), '')
        if g.has_node(left_kmer) and g.has_node(right_kmer):
            if not g.has_edge(left_kmer, right_kmer, key=gs1.genome_id):
                st1 = ('e', 'w')
                g.add_edge(left_kmer, right_kmer, key=gs1.genome_id, _support=0,
                           _type="ref_adj", _genotype= '', _dir = st1)
           
    return g, db_to_cl


def conn_shared_segments(graph, db_to_cl,genomic_segments, adj_segments):
    seg_pos = defaultdict(list)
    for n in graph.nodes:
        pos = graph.nodes[n]['_coordinate_tuple']
        seg_pos[(pos[0],pos[1], -1)].append(n)
        seg_pos[(pos[0],pos[2],1)].append(n)
    
    for pos, nodes in seg_pos.items():
        if len(nodes) == 1:
            continue
        dbls = []
        nodes = list(set(nodes))
        for n in nodes:
            dbls += db_to_cl[n]
        dbls = list(set(dbls))
        pos_ls = []
        pos_ls2 = []
        supp1 = 0
        supp2 = []
        for db in dbls:
            if db.bp_1.position == pos[1] and db.bp_1.ref_id == pos[0]:
                pos_ls.append(db)
                supp1 += db.supp
            elif db.bp_2.position == pos[1] and db.bp_2.ref_id == pos[0]:
                pos_ls.append(db)
                supp1 += db.supp
            else:
                pos_ls2.append(db)
                supp2.append(db.supp)
        if not supp1 or not supp2:
            continue
        gslist = []
        if abs(supp1 - sum(supp2)) < min(supp1, sum(supp2)):
            for db in set(pos_ls + pos_ls2):
                gslist += [gs for gs in genomic_segments[db] if pos[1] in (gs.pos1, gs.pos2)]
            for (a,b) in zip(gslist[:-1], gslist[1:]):
                adj_segments.append((a,b))
 

def _node_to_str(node_dict):
    
    if node_dict["_phase_switch"]:
        return node_dict["_coordinate"]
    
    k_THR = 2000
    lab = node_dict["_coordinate"]
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
        for subgr_num, (_,_,_type,_, cc,_) in enumerate(components_list):
            if _type == 'simple' or _type == 'indel':
                continue
            
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

def output_clusters_csv(db_to_cl, connected_components, out_file):
    db_list = defaultdict(list) 
    for subgr_num, (_,_,_type,_,cc,_) in enumerate(connected_components):
        if _type == 'simple':
            continue
        if _type == 'indel':
            for db in cc:
                db.cluster_id = "severus_" + str(subgr_num)
                db_list[subgr_num].append(db)
        else:  
            for node_id in cc:
                for db in db_to_cl[node_id]:
                    db.cluster_id = "severus_" + str(subgr_num)
                    db_list[subgr_num].append(db)
                
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\tadj_1-tadj_2\tread_support\tvaf\tgenotype\tsv_type\tdet_sv_type\tdirection\tbnd\tphaseblock\tgenome_ids\n")            
        
        for db_ls in db_list.values():
            db_ls = list(set(db_ls))
            clusters = defaultdict(list) 
            
            for br in db_ls:
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
                det_sv_type = cl[0].sv_type
                sv_type = cl[0].vcf_sv_type
                dir1 = 'TH'
                if db.direction_1 == db.direction_2:
                    dir1 = 'HH/TT'
                elif db.direction_1 == 1:
                    dir1 = 'HT'
                phase1 = ''
                if db.genotype == 'het':
                    phase1 = db.phaseset_id
                chr_type = 'intra' if db.bp_1.ref_id == db.bp_2.ref_id else 'inter'
                fout.write(f"{db.cluster_id}\t{label_1}\t{read_support}\t{db.vaf}\t{genotypes}\t{sv_type}\t{det_sv_type}\t{dir1}\t{chr_type}\t{phase1}\t{keys_text}\n")

def output_clusters_info(connected_components, out_file):
    with open(out_file, "w") as fout:
        fout.write("#cluster_id\ttype\tJunction_type\tSVs\tSV_count\tgenome_ids\t\n") 
        for subgr_num, (sv_events,junctions,_type,_,_,genome) in enumerate(connected_components):
            if _type == 'simple':
                continue
            nsv = sum([k for k in sv_events.values()])
            svs = ','.join([ f'{key}:{val} ' for key, val in sv_events.items() if val > 0])
            junc = ','.join([ f'{key}:{val} ' for key, val in junctions.items() if val > 0])
            genome_ids = ','.join(genome)
            fout.write(f"severus_{subgr_num}\t{_type}\t{junc}\t{svs}\t{nsv}\t{genome_ids}\t\n")

def cluster_adjacencies(graph, db_to_cl, components_list, target_genomes, control_genomes):
    INDEL_THR = 1.5
    INDL_n = 5
    j_type = {'11':"HH",'-1-1': "TT",'-11': "TH",'1-1': "HT", '0-0':'Interchr'}
    
    connected = defaultdict(list)
    if not components_list:
        components_list = []
    
    for node, cl in db_to_cl.items():
        if cl[0].gr_id:
            connected[cl[0].gr_id].append(node)
    
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
        junction_type = defaultdict(int)
        db_ls = defaultdict(list)
        for node in cc:
            if graph.nodes[node]['_phase_switch']:
                continue
            db_ls[db_to_cl[node][0]].append(node)
        db_ls = list(db_ls.keys())
        clusters = defaultdict(list)
        for br in db_ls:
            clusters[br.to_string()].append(br)
            
        genome = list(set([db.genome_id for db in db_ls]))
        if len(clusters) == 1:
            db = db_ls[0]
            _cp = db.vcf_sv_type if not db.sv_type else db.sv_type
            sv_type[_cp] +=1
            _jty = '0-0' if not db.bp_1.ref_id == db.bp_2.ref_id else str(db.direction_1) + str(db.direction_2)
            junction_type[j_type[_jty]]+=1
            _type = 'simple'
            if not db.bp_1.ref_id == db.bp_2.ref_id:
                _type = 'tra'
                rank = 2
            elif db.sv_type:
                rank = 3
            else:
                rank = 0
        if len(clusters) > 1:
            _type = 'complex'
            rank = len(db_ls) + 2
            for dbs in clusters.values():
                db = dbs[0]
                _cp = db.vcf_sv_type if not db.sv_type else db.sv_type
                sv_type[_cp] +=1
                _jty = '0-0' if not db.bp_1.ref_id == db.bp_2.ref_id else str(db.direction_1) + str(db.direction_2)
                junction_type[j_type[_jty]]+=1
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
                if np.median(diff) > np.median(sv_len[seq]) * INDEL_THR or len(sv_len[seq]) < INDL_n:
                    _type = 'simple'
                    rank = 4
            components_list.append((sv_type, junction_type, _type, rank, cc, genome))
            rank_ls[rank] +=1
    
    components_list.sort(key=lambda p: p[3], reverse=True)
    
    return components_list

def html_plot(graph, adj_clusters, db_to_cl, out_dir):
    db_keys = list(db_to_cl.keys())
    DODGE = 0.05
    AXIS_OFF = 500000
    k_THR = 2000
    #cov_list = [graph.nodes[n]['_coverage'] for n in graph.nodes if graph.nodes[n]['_coverage']]
    #col_segments = ['#bfd4b8', '#7fa970', '#405538']
    for subgr_num, (_,_,_type,_,cc,_) in enumerate(adj_clusters):
        if not _type == 'complex':
            continue
        
        db_list = defaultdict(list)
        segment_list = defaultdict(list)
        for c in cc:
            if c in db_keys:
                for db in db_to_cl[c]:
                    db_list[db].append(c)
            if not graph.nodes[c]['_phase_switch']:
                pos = graph.nodes[c]['_coordinate_tuple']
                segment_list[pos[0]].append([pos[1], pos[2], c, graph.nodes[c]['_haplotype']]) 
                
        segment_list = dict(sorted(segment_list.items()))
        x = []
        y = []
        lab_list = []
        y_pos = defaultdict(list)
        cols = []
        for i,(seq, pos_list) in enumerate(segment_list.items()):
            pos_list.sort(key=lambda b:b[1])
            x1 = []
            y1 = []
            for seg in pos_list:
                if not x1 or seg[0] > x1[-2]:
                    y1 += [i,i,i, None]
                else:
                    y1 += [y1[-2] + DODGE, y1[-2] + DODGE, y1[-2] + DODGE, None]
                for db in db_to_cl[seg[2]]:
                    y_pos[db].append(((seg[0], seg[1]), y1[-2]))
                x1 += [seg[0], np.mean([seg[0],seg[1]]), seg[1], None]
                len_g = graph.nodes[seg[2]]['_length']
                if len_g > k_THR:
                    length_1k = len_g / 1000
                    len_label = f"{length_1k:.1f}kb"
                else:
                    len_label = f"L:{len_g}bp" 
                cov = graph.nodes[seg[2]]['_coverage']
                
                
                hp = ','.join(list(set(['|'.join([str(s1[0]), str(s1[1])]) for s1 in seg[3]])))
                lab = graph.nodes[seg[2]]['_coordinate'] + '<br>Coverage:' + str(cov) + '<br>Length:' + len_label + '<br>Haplotype:' + hp
                lab_list += [lab, lab, lab, None]
                
            x += x1
            y += y1
            
        colors = {'11':"#256676",'-1-1': "#a20655",'-11': "#4ea6dc",'1-1': "#f19724", '0-0':'#cdcc50'}
        x_nan = [xx for xx in x if xx]
        x_limit = [min(x_nan) - AXIS_OFF, max(x_nan) + AXIS_OFF]
        
        fig = go.Figure()
        add_legend(fig, colors)
        add_dbs(fig, db_list, y_pos, colors)
        add_segments(fig, x,y, lab_list)
        plots_layout_settings(fig, list(segment_list.keys()), x_limit, subgr_num)
        fig.write_html(out_dir + '/plots/severus_' + str(subgr_num) + ".html")

def add_legend(fig, colors):
    fig.add_trace(go.Scatter(x=[-1], y=[1], legendgroup="HH", mode = 'lines',yaxis="y5",  
                             line = dict(shape = 'spline', color = colors['11'], width= 7, dash = 'solid'),
                             name="HH"))
    fig.add_trace(go.Scatter(x=[-1], y=[1], legendgroup="TT", mode = 'lines',  yaxis="y5",
                             line = dict(shape = 'spline', color = colors['-1-1'], width= 7, dash = 'solid'),
                             name="TT"))
    fig.add_trace(go.Scatter(x=[-1], y=[1], legendgroup="TH", mode = 'lines', yaxis="y5", 
                             line = dict(shape = 'spline', color = colors['-11'], width= 7, dash = 'solid'),
                             name="TH"))
    fig.add_trace(go.Scatter(x=[-1], y=[1], legendgroup="HT", mode = 'lines', yaxis="y5", 
                             line = dict(shape = 'spline', color = colors['1-1'], width= 7, dash = 'solid'),
                             name="HT"))
    fig.add_trace(go.Scatter(x=[-1], y=[1], legendgroup="Interchr", mode = 'lines', yaxis="y5", 
                             line = dict(shape = 'spline', color = colors['0-0'], width= 7, dash = 'solid'),
                             name="Interchr"))
    

def add_segments(fig, x,y, hoverdata):
    fig.add_trace(go.Scatter(
    x=x,
    y=y,
    name='Segments',
    yaxis="y5",
    line = dict(shape = 'spline', color = '#7fa970', width= 15, dash = 'solid'),
    mode='lines',
    opacity=0.9,
    showlegend=False,
    text=hoverdata,
    hoverinfo="text"))
    
    
def add_dbs(fig, db_list, y_pos, colors):
    DODGE = 0.02
    db_ls = []
    for db, nodes in db_list.items():
        if db in db_ls:
            continue
        dbs = [db2 for db2, nodes2 in db_list.items() if nodes == nodes2]
        db_ls += dbs
        dbs = list(set(dbs))
        db = dbs[0]
        
        gen_ids = ','.join(list(set([db.genome_id for db in dbs])))
        haplotypes = ','.join(list(set([str(db.haplotype_1) + '|' + str(db.haplotype_2) for db in dbs])))
        supp = ','.join([str(db.supp) for db in dbs])
        
        lab = ['supp:' + supp + '<br>sample:' + gen_ids + '<br>haplotypes:' + haplotypes + '<br>' + db.to_string()]
        
        c1 = '0-0'
        if db.bp_1.ref_id == db.bp_2.ref_id:
            c1 = str(db.direction_1) + str(db.direction_2)
            
        col = colors[c1]
        x_b = [db.bp_1.position, np.mean([db.bp_1.position, db.bp_2.position]),db.bp_2.position]
        
        pos_list = y_pos[db]
        y0 = [b for (a,b) in pos_list if db.bp_1.position in a][0]
        y2 = [b for (a,b) in pos_list if db.bp_2.position in a][0]
        if y0 == y2:
            y1 = y0 + (DODGE * 2)
        else:
            y1 = np.quantile([y0,y2], 0.75)
        y_b = [y0, y1, y2]
        
        fig.add_trace(go.Scatter(
        x=x_b,
        y=y_b,
        text=[''],
        yaxis="y5",
        line = dict(shape = 'spline', color = col, width= 2, dash = 'solid'),
        mode='lines',
        opacity=0.9,
        showlegend=False,  
        hoverinfo="text"))
        
        fig.add_trace(go.Scatter(
        x=[x_b[1]],
        y=[y_b[1]],
        name='db',
        text=lab,
        yaxis="y5",
        marker=dict(size=8, symbol="diamond-wide", color=col),
        mode='markers',
        opacity=0.9,
        showlegend=False,  
        hoverinfo="text"))
    

def plots_layout_settings(fig, chr_list, x_limit, cluster_ind):
    y_limit = [-1,len(chr_list)]
    fig.update_layout(
        xaxis=dict(
            type="linear",
            showline=True,
            zeroline=True,
            linecolor = "dimgray",
            range=x_limit,
            tickfont={"color": "black", 'size':15}
            
        ),
        yaxis5=dict(
            linecolor="dimgray",
            tickmode = 'array',
            range = y_limit,
            tickvals = list(range(len(chr_list))),
            ticktext = chr_list,
            side="left",
            tickfont={"color": "black", 'size':15},
            ticks="outside",
            title="",
            titlefont={"color": "dimgray"},
            type="linear",
            showline=True,
            zeroline=True,
        ))
    
    fig.update_layout(
        template="plotly_white",
        font_family="Helvetica"
    )
    
    fig.update_layout(legend=dict(
        orientation = 'h', xanchor = "center", x = 0.45, y= 1.2))
    
    fig.update_layout(margin=dict(l=5, r=5, b=5, pad=1))
    fig.update_xaxes(tick0=0.0, rangemode="nonnegative")
    fig.update_layout(legend={'itemsizing': 'constant'})
    fig.update_layout(font_family= "Helvetica")
    
    fig.update_layout(
        title={
            'text': 'subcluster - ' + str(cluster_ind),
            'y':0.9,
            'x':0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family = "Helvetica",
        font_color = "black",
        font_size = 15,
        title_font_family = "Helvetica",
        title_font_color = "black",
        legend_font_size = 15
    )
    
    if len(chr_list) < 3:
        height = 400
        width = 800
    elif len(chr_list) < 7:
        height = 600
        width = 1200
    else:
        height = 800
        width = 1200
    
    xlen = x_limit[1]-x_limit[0]
    if xlen > 100000000:
        width += 200
    elif xlen > 150000000:
        width +=400
     
    fig.update_layout(
        width=width,
        height=height,
       )

    
def build_breakpoint_graph(genomic_segments, adj_segments,
                           components_list, target_genomes, control_genomes):
    graph, db_to_cl = build_graph(genomic_segments, adj_segments)
    conn_shared_segments(graph, db_to_cl,genomic_segments, adj_segments)
    graph, db_to_cl = build_graph(genomic_segments, adj_segments)
    adj_clusters = cluster_adjacencies(graph, db_to_cl, components_list, target_genomes, control_genomes)
    return graph, adj_clusters, db_to_cl
            
            
def output_graphs(db_list, coverage_histograms, thread_pool, target_genomes, control_genomes, ref_lengths, args):
    keys = ['germline', 'somatic']
    
    for key in keys:
        if not key in list(db_list.keys()):
            continue
        double_breaks = db_list[key]
        
        sub_fol = "all_SVs" if key == 'germline' else 'somatic_SVs'
        out_folder = os.path.join(args.out_dir, sub_fol)
        if not os.path.isdir(out_folder):
            os.mkdir(out_folder)
            
        out_folderhtml = os.path.join(out_folder , 'plots')
        if not os.path.isdir(out_folderhtml):
            os.mkdir(out_folderhtml)
            
        all_ids = target_genomes + control_genomes if key == 'germline' else target_genomes   
            
        logger.info(f"Preparing outputs for {sub_fol}")
        
        out_breakpoint_graph = os.path.join(out_folder, "breakpoint_graph.gv")
        out_clustered_breakpoints = os.path.join(out_folder, "breakpoint_clusters.tsv")
        out_cluster_list = os.path.join(out_folder, "breakpoint_clusters_list.tsv")
        
        logger.info("\tComputing segment coverage")
        
        for db in double_breaks:
            db.cluster_id = 0
        
        (genomic_segments, adj_segments) = get_genomic_segments(double_breaks, coverage_histograms, args.phase_vcf,
                                                                key, ref_lengths, args.min_ref_flank, args.max_genomic_len, args.min_sv_size)
        components_list = []
        if key == 'germline':
            components_list = cluster_indels(double_breaks)
        logger.info("\tPreparing graph")
        graph, adj_clusters, db_to_cl = build_breakpoint_graph(genomic_segments, adj_segments, components_list, target_genomes, control_genomes)
        html_plot(graph, adj_clusters, db_to_cl, out_folder)
        #output_clusters_graphvis(graph, adj_clusters, out_breakpoint_graph)
        output_clusters_csv(db_to_cl, adj_clusters, out_clustered_breakpoints)
        
        output_clusters_info(adj_clusters, out_cluster_list)
        
        logger.info("\tWriting vcf")
        write_to_vcf(double_breaks, all_ids, out_folder, key, ref_lengths, args.no_ins)
        
            
        
    
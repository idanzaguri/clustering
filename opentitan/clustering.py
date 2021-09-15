import pickle
from igraph import *
import os.path
from os import path

import dendrogram_utils
import igraph_utils
import hal_utils

from dendrogram_utils import *
from igraph_utils import *
from hal_utils import *

igraph_utils.netlist = netlist
hal_utils.netlist = netlist

#from cophenetic_module import *
from clustering_module import *

import numpy as np
import scipy.cluster.hierarchy
from scipy.sparse import csr_matrix
from sknetwork.hierarchy import *


GEPHI_EXPORT  = False # export to Gephi
FFDG          = False # use flip flop dependency graph instead
USE_ORIG_HIER = False

# OpenTitan
if netlist.top_module.type == "top_earlgrey_KmacEnMasking1":
    level_0_hier_list = ['u_aes', 'u_alert_handler', 'u_clkmgr', 'u_csrng', 'u_dm_top', 'u_edn.', 'u_entropy_src', 'u_flash_ctrl', 'u_flash_eflash', 'u_gpio', 'u_hmac', 'u_keymgr', 'u_kmac', 'u_lc_ctrl', 'u_nmi_gen', 'u_otbn', 'u_otp_ctrl', 'u_padctrl', 'u_pinmux', 'u_pwrmgr', 'u_ram1p.*', 'u_rstmgr', 'u_rv_core_ibex', 'u_rv_plic', 'u_rv_timer', 'u_sensor_ctrl', 'u_spi_device', 'u_sram_ctrl_main', 'u_sram_ctrl_ret', 'u_tl_adapter.*', 'u_uart.*', 'u_usbdev', 'u_xbar_main', 'u_xbar_peri']
    clusters_span = range(4,40,1)
    flat_clustering_alg         = []#'multilevel', 'spinglass', 'leading_eigenvector', 'label_propagation', 'infomap']
    hierarchical_clustering_alg = ['#fastgreedy', '#walktrap', 'Paris', 'Ward', 'LouvainHierarchy']

# PicoRV32
if netlist.top_module.type == "picorv32":
    level_0_hier_list = ['cpuregs' , 'genblk1.pcpi_mul' , 'genblk2.pcpi_div']
    clusters_span = range(2,30)
    compare_method = ['vi', 'nmi', 'split-join', 'rand', 'adjusted_rand']
    flat_clustering_alg         = ['multilevel', 'spinglass', 'leading_eigenvector', 'label_propagation', 'infomap', 'fastgreedy', 'walktrap']
    hierarchical_clustering_alg = ['fastgreedy', 'walktrap', 'Paris', 'Ward', 'LouvainHierarchy']

# SweRV EH2
if netlist.top_module.type == "eh2_swerv":
    level_0_hier_list = ['dbg', 'dec', 'dma_ctrl', 'exu', 'ifu', 'lsu', 'lsu', 'pic_ctrl_inst']
    level_0_hier_list = ['lsu', 'lsu/dccm_ctl', 'ifu', 'ifu/mem_ctl', 'pic_ctrl_inst', 'dbg', 'dma_ctrl', 'dec', 'exu', 'exu/div_e1', 'exu/i0_alu_e1', 'exu/i0_alu_e4', 'exu/i1_alu_e1', 'exu/i1_alu_e4', 'exu/mul_e1']
    clusters_span = range(2,30)
    flat_clustering_alg         = ['multilevel', 'spinglass', 'leading_eigenvector', 'label_propagation', 'infomap']
    hierarchical_clustering_alg = ['-fastgreedy', 'walktrap', '-Paris', '-Ward', '-LouvainHierarchy']
 


#
# Convert netlist to graph
# load rom file if already exists
#

if path.exists("bare_netlist_graph.gml"):
    mygraph = load("bare_netlist_graph.gml")
else:
    print("Convert netlist to igraph...")
    top = netlist.get_module_by_id(1)
    mygraph = netlist2igraph(directed=True, multiple_edges=False, inputs_as_verticies=True, top=top , clean_graph=True, debug=False)
    mygraph.save("bare_netlist_graph.gml")
    mygraph.save("bare_netlist_graph.graphml")



# find N most driving cells
N = 2
degree = [mygraph.degree(vid,OUT) for vid in mygraph.vs]
largest_vertices = sorted(range(len(degree)), key = lambda sub: degree[sub])[-N:] 
print("\n{} most driving cells:".format(N))
for vid in largest_vertices:
    print("   {}: {}".format(mygraph.vs[vid]["label"], degree[vid] ) )



print("\nOriginal netlist hierarchy:")
orig_hier = parse_netlist_clusters(print_hier=True, gl_clustering=None)
orig_gl_clustering = convert_hier_list_to_gl_clustering( list(orig_hier.keys()) )

if USE_ORIG_HIER:
    golden_clustering = make_golden_clustering(mygraph, orig_gl_clustering)

else:
    print("\nLevel 0 hierarchy:")
    expected_gl_clustering = convert_hier_list_to_gl_clustering(level_0_hier_list)
    parse_netlist_clusters(print_hier=True, gl_clustering=expected_gl_clustering)
    golden_clustering = make_golden_clustering(mygraph, expected_gl_clustering)



#saving some data before converting graph to ffdg
hdl_hier_list = mygraph.vs["label"];
save_clustering(hdl_hier_list, "hdl_hier_list")

if FFDG:
    print(" ______ ______ _____   _____ ")
    print("|  ____|  ____|  __ \ / ____|")
    print("| |__  | |__  | |  | | |  __ ")
    print("|  __| |  __| | |  | | | |_ |")
    print("| |    | |    | |__| | |__| |")
    print("|_|    |_|    |_____/ \_____|")                             
    mygraph = convert2ffdg(mygraph)

mygraph = mygraph.as_undirected()




"""
t = HDLTree()
for vid,cell in enumerate(hdl_hier_list):
    if mygraph.vs[vid]['gateid'] == -1:
        break;
    t.add_cell(cell)
t.print()
"""

"""
top = netlist.get_top_module()
gates = top.get_gates(lambda g: not (g.is_vcc_gate() or g.is_gnd_gate()) )
num_of_gates = len(gates)
submodules = top.get_submodules()
num_of_submodules = len(submodules)
inputs_nets = top.get_input_nets()
for n in top.get_input_nets():
    if n.get_source() and (n.get_source().gate.is_gnd_gate() or n.get_source().gate.is_vcc_gate()):
        inputs_nets.remove(n)
    num_of_inputs = len(inputs_nets)
print("num_of_gates="+str(num_of_gates))
print("num_of_submodules="+str(num_of_submodules))
print("num_of_inputs="+str(num_of_inputs))   
"""





#  _   _ _                         _     _           _   _____ _           _            _             
# | | | (_)                       | |   (_)         | | /  __ \ |         | |          (_)            
# | |_| |_  ___ _ __ __ _ _ __ ___| |__  _  ___ __ _| | | /  \/ |_   _ ___| |_ ___ _ __ _ _ __   __ _ 
# |  _  | |/ _ \ '__/ _` | '__/ __| '_ \| |/ __/ _` | | | |   | | | | / __| __/ _ \ '__| | '_ \ / _` |
# | | | | |  __/ | | (_| | | | (__| | | | | (_| (_| | | | \__/\ | |_| \__ \ ||  __/ |  | | | | | (_| |
# \_| |_/_|\___|_|  \__,_|_|  \___|_| |_|_|\___\__,_|_|  \____/_|\__,_|___/\__\___|_|  |_|_| |_|\__, |
#                                                                                                __/ |
#                                                                                               |___/ 
print("\nComparing different hierarchical clustering algorithms using cophenetic correlation metric:")


stop_split = [2000,1000,500]
exp_hier   = [0,2,8,16]

adjacency = mygraph.get_adjacency_sparse() # for Network package
print(type(adjacency))
save_clustering(adjacency, "adjacency")



for alg in hierarchical_clustering_alg:
    if alg == 'fastgreedy':
        dendrogram = mygraph.community_fastgreedy(weights=None)
    elif alg == 'walktrap':
        dendrogram = mygraph.community_walktrap(weights=None, steps=4)
    elif alg == 'Paris':
        dendrogram = Paris().fit_transform(adjacency)
    elif alg == 'Ward':
        dendrogram = Ward().fit_transform(adjacency)
    elif alg == 'LouvainHierarchy':
        dendrogram = LouvainHierarchy().fit_transform(adjacency)
    else:
        continue
    
    print(alg)

    save_clustering(dendrogram, "dendrogram_" + alg)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    continue # FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if isinstance(dendrogram, VertexDendrogram):
        cpp_dendrogram = dendrogram.merges
    else:
        cpp_dendrogram = network_dendrogram_2_igraph_dendrogram(dendrogram) # TODO: move it to cpp code

    #for h in exp_hier:
        #for s in stop_split:
            #set_1(s)
            #set_2(h)
    print("cophenetic_correlation: (%d,%d) = %f "%(h,s,cophenetic_correlation(cpp_dendrogram, hdl_hier_list, 4, 0)))

    #print(mygraph.community_fastgreedy(weights=None).as_clustering(3).membership)
    #print(type(mygraph.community_fastgreedy(weights=None).as_clustering(3).membership))
    
    # if algorithm is also in flat clustering list, check result as flat clustering
    if alg in flat_clustering_alg:
        list_of_clustering = {}
        for method in compare_method:   
            cmp_res = [compare_communities(golden_clustering, dendrogram.as_clustering(n), method) for n in clusters_span]
            print("%-15s (%2d) %.3f"%(method, cmp_res.index(max(cmp_res))+2, max(cmp_res) ))
            list_of_clustering[cmp_res.index(max(cmp_res))+2] = 1
        if GEPHI_EXPORT:
            for clustering_size in list_of_clustering.keys():
                export_compressed_graph_to_gephi(mygraph, dendrogram.as_clustering(clustering_size), expected_gl_clustering, alg + str(clustering_size))
                export_to_gephi                 (mygraph, dendrogram.as_clustering(clustering_size), expected_gl_clustering, alg + str(clustering_size))
       

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
assert False," stop here. no need to chcek flat clustering for now" # FIXME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# ______ _       _     _____ _           _            _             
# |  ___| |     | |   /  __ \ |         | |          (_)            
# | |_  | | __ _| |_  | /  \/ |_   _ ___| |_ ___ _ __ _ _ __   __ _ 
# |  _| | |/ _` | __| | |   | | | | / __| __/ _ \ '__| | '_ \ / _` |
# | |   | | (_| | |_  | \__/\ | |_| \__ \ ||  __/ |  | | | | | (_| |
# \_|   |_|\__,_|\__|  \____/_|\__,_|___/\__\___|_|  |_|_| |_|\__, |
#                                                              __/ |
#                                                             |___/
print("\nComparing different flat clustering algorithms using different metrics:")

if 'leading_eigenvector' in flat_clustering_alg:
    print("community_leading_eigenvector")
    clustering = [mygraph.community_leading_eigenvector(n, weights=None) for n in clusters_span]
    save_clustering(clustering, "leading_eigenvector") 
    list_of_clustering = {}
    for method in compare_method:    
        cmp_res = [compare_communities(golden_clustering, c, method) for c in clustering]
        print("%-15s (%2d) %.3f"%(method, clusters_span[cmp_res.index(max(cmp_res))], max(cmp_res) ))
        list_of_clustering[ clusters_span[cmp_res.index(max(cmp_res))] ] = cmp_res.index(max(cmp_res))
    if GEPHI_EXPORT:
        for clustering_size in list_of_clustering.keys():
            export_compressed_graph_to_gephi(mygraph, clustering[list_of_clustering[clustering_size]], expected_gl_clustering, "leading_eigenvector_" + str(clustering_size))
            export_to_gephi                 (mygraph, clustering[list_of_clustering[clustering_size]], expected_gl_clustering, "leading_eigenvector_" + str(clustering_size))
     


if 'multilevel' in flat_clustering_alg:     
    print("community_multilevel")
    clustering = mygraph.community_multilevel(weights=None, return_levels=True)
    clustering = [c for c in clustering if len(c) <= clusters_span[-1]*2 ] # filter levels with large number of clusters
    save_clustering(clustering, "multilevel") 
    for method in compare_method:
        cmp_res = [compare_communities(golden_clustering, c, method) for c in clustering]
        c_index = cmp_res.index(max(cmp_res))
        num_of_clusters = len( clustering[c_index] )
        print("%-15s (%2d) %.3f"%(method, num_of_clusters, max(cmp_res) ))
        if GEPHI_EXPORT:
            export_compressed_graph_to_gephi(mygraph, clustering[c_index], expected_gl_clustering, "multilevel_" + str(num_of_clusters))
            export_to_gephi                 (mygraph, clustering[c_index], expected_gl_clustering, "multilevel_" + str(num_of_clusters))
    

for alg in flat_clustering_alg:
    if alg == 'spinglass':
        clustering = mygraph.community_spinglass(weights=None, spins=25, parupdate=False, start_temp=1, stop_temp=0.01, cool_fact=0.99, update_rule="config", gamma=1, implementation="orig", lambda_=1)
    elif alg == 'label_propagation':
        clustering = mygraph.community_label_propagation(weights=None, initial=None, fixed=None)
    elif alg == 'infomap':
        clustering = mygraph.community_infomap(edge_weights=None, vertex_weights=None, trials=10)
    else:
        continue

    print(alg)

    save_clustering(clustering, alg)

    for method in compare_method:
        cmp_res = compare_communities(golden_clustering, clustering, method) 
        print("%-15s (%2d) %.3f"%(method, len(clustering), cmp_res ))

    if GEPHI_EXPORT:
        export_compressed_graph_to_gephi(mygraph, clustering, expected_gl_clustering, alg + "_" + str(len(clustering)) )
        export_to_gephi                 (mygraph, clustering, expected_gl_clustering, alg + "_" + str(len(clustering)) )





   







"""
assert False, "----------- END ------- no need to run FFDG"

## FFDG
print("\nComparing different clustering algorithms using different metrics on FFD Graph:")

if 'fastgreedy' in clustering_alg:
    print("community_fastgreedy")
    clustering = ffdg.community_fastgreedy(weights=None)
    list_of_clustering = {}
    for method in compare_method:    
        cmp_res = [compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, clustering.as_clustering(n)), method) for n in clusters_span]
        print("%-15s (%2d) %.3f"%(method, cmp_res.index(max(cmp_res))+2, max(cmp_res) ))
        list_of_clustering[cmp_res.index(max(cmp_res))+2] = 1
   
if 'multilevel' in clustering_alg:     
    print("community_multilevel")
    clustering = ffdg.community_multilevel(weights=None, return_levels=True)
    clustering = [c for c in clustering if len(c) <= clusters_span[-1]*2 ] # filter levels with large number of clusters

    for method in compare_method:
        cmp_res = [compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, c), method) for c in clustering]
        c_index = cmp_res.index(max(cmp_res))
        num_of_clusters = len( clustering[c_index] )
        print("%-15s (%2d) %.3f"%(method, num_of_clusters, max(cmp_res) ))
     
if 'spinglass' in clustering_alg:    
    print("community_spinglass")
    clustering = ffdg.community_spinglass(weights=None, spins=25, parupdate=False, start_temp=1, stop_temp=0.01, cool_fact=0.99, update_rule="config", gamma=1, implementation="orig", lambda_=1)
    for method in compare_method:
        cmp_res = compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, clustering), method) 
        print("%-15s (%2d) %.3f"%(method, len(clustering), cmp_res ))
    
if 'walktrap' in clustering_alg:    
    print("community_walktrap")
    clustering = ffdg.community_walktrap(weights=None, steps=4)
    list_of_clustering = {}
    for method in compare_method:
        cmp_res = [compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, clustering.as_clustering(n)), method) for n in clusters_span]
        print("%-15s (%2d) %.3f"%(method, cmp_res.index(max(cmp_res))+2, max(cmp_res) ))
        list_of_clustering[cmp_res.index(max(cmp_res))+2] = 1
     
    
if 'leading_eigenvector' in clustering_alg:
    print("community_leading_eigenvector")
    clustering = [ffdg.community_leading_eigenvector(n, weights=None) for n in clusters_span]
    list_of_clustering = {}

    for method in compare_method:    
        cmp_res = [compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, c), method) for c in clustering]
        print("%-15s (%2d) %.3f"%(method, clusters_span[cmp_res.index(max(cmp_res))], max(cmp_res) ))
        list_of_clustering[ clusters_span[cmp_res.index(max(cmp_res))] ] = cmp_res.index(max(cmp_res))

if 'label_propagation' in clustering_alg:
    print("community_label_propagation")
    clustering = ffdg.community_label_propagation(weights=None, initial=None, fixed=None)

    for method in compare_method:    
        cmp_res = compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, clustering), method)
        print("%-15s (%2d) %.3f"%(method, len(clustering), cmp_res ))
 

if 'infomap' in clustering_alg:    
    print("community_infomap")
    clustering = ffdg.community_infomap(edge_weights=None, vertex_weights=None, trials=10)
    for method in compare_method:
        cmp_res = compare_communities(golden_clustering, extract_ffdg_clustering(ffdg, clustering), method) 
        print("%-15s (%2d) %.3f"%(method, len(clustering), cmp_res ))
    

"""

#split_module(mygraph, clustering, top)
#x = mygraph.subgraph(vertices=clustering[0], implementation="create_from_scratch")

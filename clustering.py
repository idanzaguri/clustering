import pickle
from igraph import *
import os.path
from os import path
import igraph_utils
import hal_utils
from igraph_utils import *
from hal_utils import *
igraph_utils.netlist = netlist
hal_utils.netlist = netlist


def save_clustering(clustering, name):
    pk_file = open(str(name) + '.pk','wb')
    pickle.dump(clustering, pk_file)
    pk_file.close()

def load_clustering(name):
    pk_file = open(str(name) + '.pk','rb')
    pk_data = pickle.load(pk_file)
    pk_file.close()
    return pk_data 

    
# OpenTitan
if netlist.top_module.type == "top_earlgrey_KmacEnMasking1":
    level_0_hier_list = ['u_aes', 'u_alert_handler', 'u_clkmgr', 'u_csrng', 'u_dm_top', 'u_edn.', 'u_entropy_src', 'u_flash_ctrl', 'u_flash_eflash', 'u_gpio', 'u_hmac', 'u_keymgr', 'u_kmac', 'u_lc_ctrl', 'u_nmi_gen', 'u_otbn', 'u_otp_ctrl', 'u_padctrl', 'u_pinmux', 'u_pwrmgr', 'u_ram1p.*', 'u_rstmgr', 'u_rv_core_ibex', 'u_rv_plic', 'u_rv_timer', 'u_sensor_ctrl', 'u_spi_device', 'u_sram_ctrl_main', 'u_sram_ctrl_ret', 'u_tl_adapter.*', 'u_uart.*', 'u_usbdev', 'u_xbar_main', 'u_xbar_peri']
    clusters_span = range(4,40,1)
    compare_method = ['vi', 'nmi', 'split-join', 'rand', 'adjusted_rand']
    clustering_alg = ['fastgreedy', 'DISABLED walktrap', 'DISABLED multilevel', 'DISABLED spinglass']

# PicoRV32
if netlist.top_module.type == "picorv32":
    level_0_hier_list = ['cpuregs' , 'genblk1.pcpi_mul' , 'genblk2.pcpi_div']
    clusters_span = range(2,30)
    compare_method = ['vi', 'nmi', 'split-join', 'rand', 'adjusted_rand']
    clustering_alg = ['fastgreedy', 'walktrap', 'multilevel', 'spinglass']

# SweRV EH2
if netlist.top_module.type == "eh2_swerv":
    level_0_hier_list = ['dbg', 'dec', 'dma_ctrl', 'exu', 'ifu', 'lsu', 'lsu', 'pic_ctrl_inst']
    level_0_hier_list = ['lsu', 'lsu/dccm_ctl', 'ifu', 'ifu/mem_ctl/icache.*', 'pic_ctrl_inst', 'dbg', 'dma_ctrl', 'dec', 'exu', 'exu/div_e1', 'exu/i0_alu_e1', 'exu/i0_alu_e4', 'exu/i1_alu_e1', 'exu/i1_alu_e4', 'exu/mul_e1']
    clusters_span = range(2,30)
    compare_method = ['vi', 'nmi', 'split-join', 'rand', 'adjusted_rand']
    clustering_alg = ['fastgreedy', 'walktrap', 'multilevel', 'spinglass']



if path.exists("bare_netlist_graph.gml"):
    mygraph = load("bare_netlist_graph.gml")
else:
    print("Convert netlist to igraph...")
    top = netlist.get_module_by_id(1)
    mygraph = netlist2igraph(directed=False, multiple_edges=False, inputs_as_verticies=True, top=top , clean_graph=True, debug=False)
    mygraph.save("bare_netlist_graph.gml")

# find N most driving cells
N = 2
degree = [mygraph.degree(vid,OUT) for vid in mygraph.vs]
largest_vertices = sorted(range(len(degree)), key = lambda sub: degree[sub])[-N:] 
print("\n{} most driving cells:".format(N))
for vid in largest_vertices:
    print("   {}: {}".format(mygraph.vs[vid]["label"], degree[vid] ) )


#print("\nOriginal netlist hierarchy:")
#orig_hier = parse_netlist_clusters(print_heir=False, gl_clustering=None)
#orig_gl_clustering = convert_hier_list_to_gl_clustering( list(orig_hier.keys()) )

print("\nLevel 0 hierarchy:")
expected_gl_clustering = convert_hier_list_to_gl_clustering(level_0_hier_list)
parse_netlist_clusters(print_heir=True, gl_clustering=expected_gl_clustering)

#golden_clustering = make_golden_clustering(mygraph, orig_gl_clustering)
golden_clustering = make_golden_clustering(mygraph, expected_gl_clustering)

print("\nComparing different clustering algorithms using different metrics:")

if 'fastgreedy' in clustering_alg:
    print("community_fastgreedy")
    clustering = mygraph.community_fastgreedy(weights=None)
    save_clustering(clustering, "fastgreedy") 
    list_of_clustering = {}
    for method in compare_method:    
        cmp_res = [compare_communities(golden_clustering, clustering.as_clustering(n), method) for n in clusters_span]
        print("%-15s (%2d) %.3f"%(method, cmp_res.index(max(cmp_res))+2, max(cmp_res) ))
        list_of_clustering[cmp_res.index(max(cmp_res))+2] = 1
    for clustering_size in list_of_clustering.keys():
        export_compressed_graph_to_gephi(mygraph, clustering.as_clustering(clustering_size), expected_gl_clustering, "fastgreedy_" + str(clustering_size))
        export_to_gephi                 (mygraph, clustering.as_clustering(clustering_size), expected_gl_clustering, "fastgreedy_" + str(clustering_size))
     
   
if 'multilevel' in clustering_alg:     
    print("community_multilevel")
    clustering = mygraph.community_multilevel(weights=None, return_levels=True)
    clustering = [c for c in clustering if len(c) <= clusters_span[-1]*2 ] # filter levels with large number of clusters
    save_clustering(clustering, "multilevel") 

    for method in compare_method:
        cmp_res = [compare_communities(golden_clustering, c, method) for c in clustering]
        c_index = cmp_res.index(max(cmp_res))
        num_of_clusters = len( clustering[c_index] )
        print("%-15s (%2d) %.3f"%(method, num_of_clusters, max(cmp_res) ))
        export_compressed_graph_to_gephi(mygraph, clustering[c_index], expected_gl_clustering, "multilevel_" + str(num_of_clusters))
        export_to_gephi                 (mygraph, clustering[c_index], expected_gl_clustering, "multilevel_" + str(num_of_clusters))
     

if 'spinglass' in clustering_alg:    
    print("community_spinglass")
    clustering = mygraph.community_spinglass(weights=None, spins=25, parupdate=False, start_temp=1, stop_temp=0.01, cool_fact=0.99, update_rule="config", gamma=1, implementation="orig", lambda_=1)
    save_clustering(clustering, "spinglass") 
    for method in compare_method:
        cmp_res = compare_communities(golden_clustering, clustering, method) 
        print("%-15s (%2d) %.3f"%(method, len(clustering), cmp_res ))
    export_compressed_graph_to_gephi(mygraph, clustering, expected_gl_clustering, "spinglass_" + str(len(clustering)) )
    export_to_gephi                 (mygraph, clustering, expected_gl_clustering, "spinglass_" + str(len(clustering)) )
    

if 'walktrap' in clustering_alg:    
    print("community_walktrap")
    clustering = mygraph.community_walktrap(weights=None, steps=4)
    save_clustering(clustering, "walktrap")     
    list_of_clustering = {}
    for method in compare_method:
        cmp_res = [compare_communities(golden_clustering, clustering.as_clustering(n), method) for n in clusters_span]
        print("%-15s (%2d) %.3f"%(method, cmp_res.index(max(cmp_res))+2, max(cmp_res) ))
        list_of_clustering[cmp_res.index(max(cmp_res))+2] = 1
    for clustering_size in list_of_clustering.keys():
        export_compressed_graph_to_gephi(mygraph, clustering.as_clustering(clustering_size), expected_gl_clustering, "walktrap_" + str(clustering_size))
        export_to_gephi                 (mygraph, clustering.as_clustering(clustering_size), expected_gl_clustering, "walktrap_" + str(clustering_size))
     

#split_module(mygraph, clustering, top)
#x = mygraph.subgraph(vertices=clustering[0], implementation="create_from_scratch")

# community_leading_eigenvector
#if clustering_alg["leading_eigenvector"]:
#    nmi = [compare_communities(golden_clustering, mygraph.community_leading_eigenvector(n, weights=None), "nmi") for n in clusters_span]
#    print(  len(mygraph.community_leading_eigenvector(56, weights=None) )  )
#    print( "community_leading_eigenvector({}) nmi={}".format(nmi.index(max(nmi))+2, max(nmi)) )
#    clustering_results["leading_eigenvector"] = {"clustering": 1, "num_of_clusters": nmi.index(max(nmi))+2}
# community_label_propagation
#if clustering_alg["label_propagation"]:
#    clustering = mygraph.community_label_propagation(weights=None, initial=None, fixed=None)
#    nmi = compare_communities(golden_clustering, clustering, "nmi") 
#    print( "community_label_propagation({}) nmi={}".format(len(clustering), nmi ) )
#    #clustering_report(mygraph, clustering, expected_gl_clustering)
#    clustering_results["label_propagation"] = clustering

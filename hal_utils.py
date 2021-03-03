from igraph import *
import re
import copy
import math

def split_module(graph, clustering, module):
    for i, set in enumerate(clustering, start=0):
        gates_id = [ graph.vs[vid]["gateid"] for vid in set ] 
        gates = [ netlist.get_gate_by_id( id ) for id in gates_id if id != -1 ]
        candidate = netlist.create_module(module.name+"_"+str(i), module, gates)


# top is represented via ''
def get_gate_cluster(gate, gl_clustering = None, debug = 0):
    m = re.split("/",gate) 
    m = m[:-1] # remove cell from path
    gate_path = '/'.join(m)

    if gl_clustering != None:
        for hier in gl_clustering:
            #re.sub('[\[\]]','@',hier) - in order to handle with hier contains regex reserved charters
            if(re.match("^" + re.sub('[\[\]]','@',hier), re.sub('[\[\]]','@',gate_path))):
                return hier
        assert False ,"TOP ('') must be included"# FIXME
        return '' # consider cell belongs to top
    return gate_path

def parse_netlist_clusters(print_heir = True, gl_clustering = None):
    if gl_clustering != None:
        hier_dict = {cluster:0 for cluster in gl_clustering}
        hier_dict[''] = 0 # adding 'top' which is the default cluster to list after sorting the hierarchies
    else:
        hier_dict = {}
        
    for gate in netlist.get_gates(lambda g: not( g.is_gnd_gate() or g.is_vcc_gate() ) ):     
        parent = get_gate_cluster(gate.name, gl_clustering) 
        if parent not in hier_dict:
            assert gl_clustering == None, "out of the blue cluster!"
            hier_dict[parent] = 1
        else:
            hier_dict[parent] += 1

    if print_heir:
        largest_path = max(hier_dict.keys(), key=len)
        padding = len(largest_path)

        print("%s | #Cells " % ( "Hierarchy".ljust(padding)))
        print("%s-+---------" % ("-"*padding) )
        num_of_cells = 0
        for path in sorted(hier_dict.keys()):
            print("%s | %7d "  % ( path.ljust(padding), hier_dict[path] ))
            num_of_cells += hier_dict[path]
            assert hier_dict[path] != 0 ,"empty hier!"
        print("%s-+---------" % ("-"*padding) )
        print("%s | %7d" % ( "Total".ljust(padding) , num_of_cells))    
            
    return hier_dict


def convert_hier_list_to_gl_clustering(hier_list):
    if '' in hier_list:
        hier_list.remove('')
    gl_clustering = sorted(hier_list, key = lambda p: len(re.split('/',p)), reverse=True) # sort by hier depth
    gl_clustering.append('')
    return gl_clustering
    

def make_golden_clustering(graph, gl_clustering):
    cluster2id = {}
    for i in range(len(gl_clustering)):
        cluster2id[ gl_clustering[i] ] = i
    return Clustering( [ cluster2id[get_gate_cluster(vid, gl_clustering)] for vid in graph.vs["label"] ] )


def clustering_report(graph, clustering, gl_clustering=None):
    for c in sorted(gl_clustering):
        print("%-3s" % c, end = ",")
    print()

    for i, cluster in enumerate(clustering):
        histogram = {key:0 for key in gl_clustering}
        for vid in cluster:
            hier = get_gate_cluster(graph.vs[vid]["label"], gl_clustering)
            histogram[hier] += 1
        #print(histogram)
        print("%3d)"%(i), end = " ")
        for bin in sorted(histogram.keys()):
            print("%4d" % histogram[bin], end = " | ")
        print()
        
            
def export_to_gephi(g, clustering, gl_clustering, filename):
    
    # This code could potentially helps some layout algorithm to converge into the selected 'clustering'.
    # However, some fine tuning is needed
    # cluster_crossing_edges = clustering.crossing()
    # graph.es["dummy_weight"] = [1+9999*cross for cross in cluster_crossing_edges]
 
    num_of_clusters = len(clustering)
    largest_cluster = max( [ len(cluster) for cluster in clustering ] )

    intra_cluster_size = math.ceil(math.sqrt(largest_cluster))
    inter_cluster_size = math.ceil(math.sqrt(num_of_clusters))
    margin = intra_cluster_size*3

    for i, cluster in enumerate(clustering):
        cluster_pos = (int(i/inter_cluster_size), i%inter_cluster_size)
        for j, vid in enumerate(cluster):
            v_pos = (int( (j)/intra_cluster_size), (j)%intra_cluster_size)
            g.vs[vid]["ypos"] = cluster_pos[0]*margin + v_pos[0]
            g.vs[vid]["xpos"] = cluster_pos[1]*margin + v_pos[1]
            g.vs[vid]["orig_cluster"] = get_gate_cluster(g.vs[vid]["label"])
            g.vs[vid]["expected_cluster"] = get_gate_cluster(g.vs[vid]["label"], gl_clustering)
            g.vs[vid]["cluster"] = i
    g.save(filename + ".gml")


def export_compressed_graph_to_gephi(graph, clustering, gl_clustering, filename):
    g = copy.deepcopy(graph)

    g.vs["expected_cluster"] = [ get_gate_cluster(g.vs[v.index]["label"], gl_clustering) for v in g.vs ]

    new_clustering = []
    new_vid = 0
    new_vertices = [-1]*g.vcount()


    counter = 0;
    total = g.vcount();

    for i, cluster in enumerate(clustering):
        counter+=len(cluster)
        print("@@@@ " + str(i) + " (+" + str(len(cluster)) + ")   " + str(counter) + "/" + str(total))

        exp_clusters_in_set = set(g.vs[v]["expected_cluster"] for v in cluster)
        #exp_clusters_in_set = set( get_gate_cluster(g.vs[v]["label"], gl_clustering) for v in cluster)

        exp_cluster_to_new_vid = {cluster:i for i,cluster in enumerate(exp_clusters_in_set, start=new_vid)}
        new_vid += len(exp_cluster_to_new_vid)
        new_clustering.append(new_vid)

        for vid in cluster:
            new_vertices[vid] = exp_cluster_to_new_vid[g.vs[vid]["expected_cluster"]]
            
    assert -1 not in new_vertices, "for some reason at least one vertex was not compressed"
        
        
    g.vs["size"] = 20
    if not g.is_weighted():
        g.es["weight"] = 1
    
                       
    #g.vs["label"] = g.vs["expected_cluster"]                                           
    g.contract_vertices(new_vertices ,combine_attrs=dict(label = "first", size = "sum"))
    g.simplify(combine_edges=dict(weight="sum "))

    num_of_clusters = len(new_clustering)
    largest_cluster = max([ (new_clustering[i+1] - new_clustering[i]) for i in range(num_of_clusters-1) ] )
    
    intra_cluster_size = math.ceil(math.sqrt(largest_cluster))
    inter_cluster_size = math.ceil(math.sqrt(num_of_clusters))
    margin = intra_cluster_size*2

    last_index = 0
    for i, index in enumerate(new_clustering):
    
        cluster_pos = (int(i/inter_cluster_size), i%inter_cluster_size)

        for vid in range(last_index, index):
            v_pos = (int( (vid-last_index)/intra_cluster_size), (vid-last_index)%intra_cluster_size)
            g.vs[vid]["ypos"] = cluster_pos[0]*margin + v_pos[0]
            g.vs[vid]["xpos"] = cluster_pos[1]*margin + v_pos[1]
            g.vs[vid]["orig_cluster"] = get_gate_cluster(g.vs[vid]["label"])
            g.vs[vid]["expected_cluster"] = get_gate_cluster(g.vs[vid]["label"], gl_clustering)
            g.vs[vid]["cluster"] = i
        last_index = index

    hier = parse_netlist_clusters(print_heir=False, gl_clustering=gl_clustering)
    for v in g.vs:
        v["size_scaled"] = v["size"] / hier[v["expected_cluster"]]
        
    g.save(filename + "_compressed.gml")

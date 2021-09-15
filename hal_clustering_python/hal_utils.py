from igraph import *
import pickle
import re
import copy
import math


def save_clustering(clustering, name):
    pk_file = open("pk_files/" + str(name) + '.pk','wb')
    pickle.dump(clustering, pk_file)
    pk_file.close()

def load_clustering(name):
    pk_file = open(str(name) + '.pk','rb')
    pk_data = pickle.load(pk_file)
    pk_file.close()
    return pk_data 

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

def parse_netlist_clusters(print_hier = True, gl_clustering = None, top = None):
    if top == None:
        top = netlist.top_module
         
    assert len(top.submodules) == 0, "Currently there is no support for modules with  submodules"
 
    if gl_clustering != None:
        hier_dict = {cluster:0 for cluster in gl_clustering}
        hier_dict[''] = 0 # adding 'top' which is the default cluster to list after sorting the hierarchies
    else:
        hier_dict = {}
        
    for gate in top.get_gates(lambda g: not( g.is_gnd_gate() or g.is_vcc_gate() ) ):     
        parent = get_gate_cluster(gate.name, gl_clustering) 
        if parent not in hier_dict:
            assert gl_clustering == None, "out of the blue cluster!"
            hier_dict[parent] = 1
        else:
            hier_dict[parent] += 1

    if print_hier:
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
    # adding coordinates for Gephi GeoLayout plugin (ypos,xpos)
 
    num_of_clusters = len(clustering)
    largest_cluster = max( [ len(cluster) for cluster in clustering ] )

    cluster_side_length = math.ceil(math.sqrt(largest_cluster))
    inter_cluster_size  = math.ceil(math.sqrt(num_of_clusters))

    for i, cluster in enumerate(clustering):
        cluster_pos = (int(i/inter_cluster_size), i%inter_cluster_size)
        for vertex_lpos, vid in enumerate(cluster):
            v_pos = ( int( vertex_lpos/cluster_side_length ), vertex_lpos%cluster_side_length )
            g.vs[vid]["ypos"] = cluster_pos[0]*2*cluster_side_length + v_pos[0]
            g.vs[vid]["xpos"] = cluster_pos[1]*2*cluster_side_length + v_pos[1]
            g.vs[vid]["orig_cluster"] = get_gate_cluster(g.vs[vid]["label"])
            g.vs[vid]["expected_cluster"] = get_gate_cluster(g.vs[vid]["label"], gl_clustering)
            g.vs[vid]["cluster_id"] = i
    g.save(filename + ".gml")

  


def remove_dummies( graph, clustering = None):                    
    dummy_vertices = [ v.index for v in graph.vs.select(gateid=-1)]
    graph.delete_vertices( dummy_vertices )
    if clustering:
        return Clustering ( [cluster_id for cluster_id,set in enumerate(clustering) for vid in set if vid not in dummy_vertices] )



def export_compressed_graph_to_gephi(graph, clustering, gl_clustering, filename):
    g = copy.deepcopy(graph)  
    clustering = remove_dummies(g, clustering)

    #
    # In each cluster, collapse all vertices with the same gl_cluster into one vertex
    #
 
    # list of sorted numbers- L[a,b,c,...] where each compressed cluster i'th contains vids in [L[i-1],[i])  {L[0-1] = 0}
    new_clustering = []
    new_vid = 0
    new_vertices = [-1]*g.vcount() # map original vertices to the compressed ones

    g.vs["expected_cluster"] = [ get_gate_cluster(g.vs[v.index]["label"], gl_clustering) for v in g.vs ]

    for i, cluster in enumerate(clustering):                
        
        # detect all expected clusters with elements in current cluser and make map from expected_cluster to new vid
        exp_clusters_in_set = set(g.vs[v]["expected_cluster"] for v in cluster) 
        exp_cluster_to_new_vid = {cluster:i for i,cluster in enumerate(exp_clusters_in_set, start=new_vid)}
        new_vid += len(exp_cluster_to_new_vid)
        new_clustering.append(new_vid)

        for vid in cluster:
            new_vertices[vid] = exp_cluster_to_new_vid[g.vs[vid]["expected_cluster"]]
            
    assert -1 not in new_vertices, "for some reason at least one vertex was not compressed"

    g.vs["size"] = 1
    if not g.is_weighted():
        g.es["weight"] = 1
                       
    g.contract_vertices(new_vertices ,combine_attrs=dict(label = "first", size = "sum"))
    g.simplify(combine_edges=dict(weight="sum"))


    # adding coordinates for Gephi GeoLayout plugin (ypos,xpos)
    num_of_clusters = len(new_clustering)
    largest_cluster = max([ (new_clustering[i+1] - new_clustering[i]) for i in range(num_of_clusters-1) ] )
    
    cluster_side_length = math.ceil(math.sqrt(largest_cluster))
    inter_cluster_size  = math.ceil(math.sqrt(num_of_clusters))

    #print(cluster_side_length) # 3
    #print(inter_cluster_size)  # 4
    
    last_index = 0
    for i, index in enumerate(new_clustering):
        cluster_pos = (int(i/inter_cluster_size), i%inter_cluster_size)

        for vid in range(last_index, index):
            vertex_lpos = vid-last_index # linear position
            v_pos = ( int( vertex_lpos/cluster_side_length ), vertex_lpos%cluster_side_length )

            g.vs[vid]["ypos"] = cluster_pos[0]*2*cluster_side_length + v_pos[0]
            g.vs[vid]["xpos"] = cluster_pos[1]*2*cluster_side_length + v_pos[1]
            g.vs[vid]["expected_cluster"] = get_gate_cluster(g.vs[vid]["label"], gl_clustering)
            g.vs[vid]["label"] = g.vs[vid]["expected_cluster"]            
            g.vs[vid]["cluster_id"] = i

        last_index = index

    
    # Generate scale sized vertices (inter-cluster and intra-cluster)
    # value is in (0,1000] since Gephi can't read integer(1) and float(0.x) at the same column
    
    hier = parse_netlist_clusters(print_hier=False, gl_clustering=gl_clustering)
    for v in g.vs:
        v["inter_scaled"] = round( (v["size"] / hier[v["expected_cluster"]])*1000, 0)
        assert v["inter_scaled"] <= 1000, "Unexpected scaling result"

    for cid in range(num_of_clusters):
        cluster_size = sum( [v["size"] for v in g.vs.select(cluster_id = cid)] )
        for v in g.vs.select(cluster_id = cid):
            v["intra_scaled"] = round( (v["size"] / cluster_size )*1000, 0)
            assert v["intra_scaled"] <= 1000, "Unexpected scaling result"

    g.save(filename + "_compressed.gml")




FLIPFLOP_TYPE_LIST =  [
'FD1', 'FD1P', 'FD1S', 'FD1SP', 'FD2', 'FD2P', 'FD2S', 'FD2SP', 'FD3', 'FD3P', 'FD3S', 'FD3SP', 'FD4', 'FD4P','FD4S', 'FD4SP',
'FJK1', 'FJK1P', 'FJK1S', 'FJK1SP', 'FJK2', 'FJK2P', 'FJK2S', 'FJK2SP', 'FJK3', 'FJK3P', 'FJK3S', 'FJK3SP',
'FDS2', 'FDS2L', 'FDS2LP', 'FDS2P', 'FT2', 'FT2P', 'FT4', 'FT4P', 'FD2TS', 'FD2TSP',
'LD1', 'LD1P', 'LD2', 'LD2P', 'LD3', 'LD3P', 'LD4', 'LD4P', 'LSR0', 'LSR0P', 'LSR2', 'RAM1', 'LD1X4', 'LD1X4P']


def ffdg_traverse(vertex, new_vid):
    if vertex["visited"] or vertex["celltype"] in FLIPFLOP_TYPE_LIST:
        return new_vid

    vertex["visited"] = True;
    vertex["ffdg_id"] = new_vid;

    for p_vid in vertex.predecessors():
        ffdg_traverse(p_vid, new_vid)
    return new_vid+1



def convert2ffdg(graph):
    g = copy.deepcopy(graph)  

    g.vs["visited"] = False
    g.vs["ffdg_id"] = -1
    for v in g.vs:
        v["orig_vid"] = str(v.index) + ','

    new_vid = 0;
    for vertex in g.vs.select(lambda vertex: vertex["celltype"] in FLIPFLOP_TYPE_LIST):
        vertex["visited"] = True;
        vertex["ffdg_id"] = new_vid;
        new_vid += 1
        for p_vid in vertex.predecessors():
            new_vid = ffdg_traverse(p_vid, new_vid)
    
    #assert all(g.vs["visited"]), "There is a least one vertex which was not visited"
    for v in g.vs.select(lambda vertex: vertex["visited"] == False):
        print("ERROR: " + v["label"] + "   (" + v["celltype"] + ")")
        v["visited"] = True
        v["ffdg_id"] = new_vid;
        new_vid += 1

    new_vertices = [ v["ffdg_id"] for v in g.vs ]
    g.contract_vertices(new_vertices ,combine_attrs=dict(orig_vid = "concat"))
    g.simplify(combine_edges=dict(weight="sum"))
    return g;

def extract_ffdg_clustering(ffdg, clustering):
    m = {}
    for cluster_id, cluster in enumerate(clustering):
        for vid in cluster:
            orig_vids = [int(v) for v in ffdg.vs[vid]["orig_vid"][:-1].split(',')]
            for v in orig_vids:
                m[v] = cluster_id;
                
    return Clustering(  [m[v] for v in range(len(m)) ] )

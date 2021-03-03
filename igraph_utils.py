from igraph import *

def netlist2igraph(directed = False, multiple_edges = False , inputs_as_verticies = True, top = None, clean_graph = True, debug = False):
    if directed:
        g = Graph().as_directed()
    else:
        g = Graph()

    if not top:
        top = netlist.get_top_module()

    gates = top.get_gates(lambda g: not (g.is_vcc_gate() or g.is_gnd_gate()) )
    num_of_gates = len(gates)

    if inputs_as_verticies:
        inputs_nets = top.get_input_nets()
        for n in top.get_input_nets():
            if n.get_source() and (n.get_source().gate.is_gnd_gate() or n.get_source().gate.is_vcc_gate()):
                inputs_nets.remove(n)
        num_of_inputs = len(inputs_nets)
    else:
        num_of_inputs = 0

    if debug:
        print(top.name)
        print("num_of_gates="+str(num_of_gates))
        print("num_of_inputs="+str(num_of_inputs))

    g.add_vertices(num_of_gates + num_of_inputs)
    vid = 0

    #
    # There are two type of vertices- gates and inputs port
    # 

    # phase 1 - dealing with vertices representing real cells
    gid_to_vid = {}
    for gate in gates:
        g.vs[vid]["label"] = gate.name
        g.vs[vid]["cell_type"] = gate.type.name
        g.vs[vid]["netid"] = -1
        g.vs[vid]["gateid"] = gate.id
        assert gate.id not in gid_to_vid, "Gate_id="+str(gate.id)+" already exists in gid_to_vid!"
        gid_to_vid[gate.id] = vid
        vid += 1

    # phase 2 - dealing with vertices of inputs nets (dummies)
    nid_to_vid = {}
    if inputs_as_verticies:
        for net in inputs_nets:
            g.vs[vid]["label"] = net.get_name()
            g.vs[vid]["netid"] = net.id
            g.vs[vid]["gateid"] = -1
            nid_to_vid[net.id] = vid
            vid += 1
    
    
    edges = []

    # convert inputs nets to edges
    if inputs_as_verticies:
        for net in inputs_nets:
            for dest in net.destinations:
                #skip if dest is in another module
                if dest.gate.module.name != top.name:
                    continue

                if debug:
                    print("INPUT:" + net.name +" -> " + dest.gate.name + "." + dest.pin)
        
                v1 = nid_to_vid[net.id]
                v2 = gid_to_vid[dest.gate.id]
                e = (v1 , v2)
                if not directed:
                    e = (min(e), max(e))
                edges.append(e)


    # convert internal nets to edges
    for net in top.get_internal_nets():
        # skipping undriven ports/nets
        if net.num_of_sources == 0:
            if debug:
                print("undriven net: " + net.name)
            continue
    
        assert net.num_of_sources == 1, "Only one source per net"

        src_g = net.sources[0].gate
        if src_g.is_gnd_gate() or src_g.is_vcc_gate():
            continue

        for dest in net.destinations:
            #skip if dest is another module
            if dest.gate.module.name != top.name:
                continue

            if debug:
                print("NET: " + net.sources[0].gate.name + "." + net.sources[0].pin + " -> " + dest.gate.name + "." + dest.pin)
        
            v1 = gid_to_vid[net.sources[0].gate.id]
            v2 = gid_to_vid[dest.gate.id]
            e = (v1 , v2)
            if not directed:
                e = (min(e), max(e))
            edges.append(e)



    if multiple_edges:
        g.add_edges(edges)
    else:
        # merge identical nets (i.e. same source and dest.) assign relative weight
        seen = {}
        dupes = []
        for e in edges:
            if e not in seen:
                seen[e] = 1      
            else:
                if seen[e] == 1:
                    dupes.append(e)
                seen[e] += 1
 
        g.add_edges(seen.keys())

        if len(dupes) > 0:
            g.es["weight"] = 1 # assign default weight 1 for all edges in case that there are duplicates nets
            for e in dupes:
                eid = g.get_eid(e[0], e[1], directed, False)
                assert eid != -1
                g.es[eid]["weight"] = seen[e];


    if clean_graph:
        clusters = g.clusters(mode=WEAK)
        if len(clusters) > 1:
            #for c in clusters:
                #for v in c:
                #    print(g.vs[v]["label"])
                #assert False,""
            clusters_size = [len(c) for c in clusters]
            print("Graph clean: %d separated sub-graphs were detected. continue with the largest one - %s"%(len(clusters) , str(clusters_size)))
            g = g.subgraph( clusters[ clusters_size.index(max(clusters_size)) ] , implementation="create_from_scratch")

    return g

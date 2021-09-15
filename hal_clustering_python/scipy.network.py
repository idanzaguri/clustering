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



from IPython.display import SVG

import numpy as np
import scipy.cluster.hierarchy
from scipy.sparse import csr_matrix
from sknetwork.hierarchy import *

from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

from sknetwork.data import karate_club, painters, movie_actor
from sknetwork.hierarchy import Paris, BiParis, cut_straight, dasgupta_score, tree_sampling_divergence
from sknetwork.visualization import svg_graph, svg_digraph, svg_bigraph, svg_dendrogram


graph = karate_club(metadata=True)
adjacency = graph.adjacency
position = graph.position
paris = Paris()
dendrogram = paris.fit_transform(adjacency)


image = svg_dendrogram(dendrogram)

SVG(image)
labels = cut_straight(dendrogram)
print(labels)

n_clusters = 4
labels, dendrogram_aggregate = cut_straight(dendrogram, n_clusters, return_dendrogram=True)
print(labels)

_, counts = np.unique(labels, return_counts=True)
image = svg_dendrogram(dendrogram_aggregate, names=counts, rotate_names=False)
SVG(image)

image = svg_graph(adjacency, position, labels=labels)
SVG(image)


assert False, "dddd"

top = netlist.get_module_by_id(1)
mygraph = netlist2igraph(directed=True, multiple_edges=False, inputs_as_verticies=True, top=top , clean_graph=True, debug=False)
mygraph = mygraph.as_undirected()



#gates = [v["label"] for v in mygraph.vs]

#print(gates)


print(mygraph.vcount())
clustering = mygraph.community_fastgreedy(weights=None)
print(type(clustering))

adjacency = mygraph.get_adjacency_sparse()
paris = Paris()

dendrogram = paris.fit_transform(adjacency)
print(dendrogram.shape)
assert False,"dddd"
print( np.round(dendrogram,2) )
for i in range(4):
    print( dendrogram[i*100][0] )
    print( dendrogram[i*100][1] )
    print( dendrogram[i*100][2] )
    print( dendrogram[i*100][3] )

#scipy.cluster.hierarchy.dendrograme(dendrogram)

"""
#X = [[i] for i in [2, 8, 0, 4, 1, 9, 9, 0]]
X = paris.fit_transform(adjacency)
Z = linkage(X, 'ward')
fig = plt.figure(figsize=(25, 10))
dn = dendrogram(Z)


Z = linkage(X, 'single')
fig = plt.figure(figsize=(25, 10))
dn = dendrogram(Z)
plt.show()
"""

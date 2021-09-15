import math
import re
import pickle
from igraph import *
import os.path
from os import path


import dendrogram_utils
from dendrogram_utils import *
from clustering_module import *


import numpy as np
import scipy.cluster.hierarchy
from scipy.sparse import csr_matrix
from sknetwork.hierarchy import *







design = "eh2_swerv"

pk_file = open('pk_files/' + design +  '/adjacency.pk','rb')
adjacency = pickle.load(pk_file)
pk_file.close()

#dendrogram = Ward().fit_transform(adjacency)
#dendrogram = Paris().fit_transform(adjacency)
#dendrogram = LouvainHierarchy().fit_transform(adjacency)

pk_file = open('pk_files/' + design +  '/Paris_dendrogram.pk','rb') #dendrogram_fastgreedy
dendrogram = pickle.load(pk_file)
pk_file.close()

pk_file = open('pk_files/' + design +  '/hdl_hier_list.pk','rb')
hdl_hier_list = pickle.load(pk_file)
pk_file.close()





# Recovering the real HDL tree using hdl_hier_list
t = HDLTree()

num_dummy_vid = len(hdl_hier_list)
for vid,cell in enumerate(hdl_hier_list):
	if re.match("^NET___", cell):
		num_dummy_vid -= vid
		break
	t.add_cell(cell, vid)

num_dummy_vid = 0
dendo = MyDendrogram(dendrogram, num_dummy_vid)

assert len(hdl_hier_list) == len(dendo.dendrogram)+1, "ERROR!"
print("NUM_ELEMENTS " + str(dendo.num_elements))



#
# this code use to check how good is the dendrogram in terms of balance
#
if False:
	merge = dendo.root;
	for i in range(10):
		split_points = dendo.split(10, merge);
		for s in split_points:
			if dendo.is_element(s):
				print("CELL! " + str(s))
			else:
				print(str(len(dendo.collapse_merge(s))) + " " + str(s))
				merge = s
		print("------------------")





stop_split = [800]#[2000,1000,500]
exp_hier   = [0,5,10] #[0,5,8,16]

if isinstance(dendrogram, VertexDendrogram):
    cpp_dendrogram = dendrogram.merges
else:
    cpp_dendrogram = network_dendrogram_2_igraph_dendrogram(dendrogram) # TODO: move it to cpp code



print("------------------------------------------")
for h in exp_hier:
    for s in stop_split:
        print("cophenetic_correlation: " + str(h) +"," + str(s))
        set_MAX_CLUSTER_SIZE(s)
        set_EXP_NUM_HIER(h)
        print("cophenetic_correlation: (%d,%d) = %f "%(h,s,cophenetic_correlation(cpp_dendrogram, hdl_hier_list, 4, 1)))

    



sys.exit()





















level_coeff = 1
def cluster_entropy(cells_paths, level = 0):
	#print(level)
	#print(cells_paths)
	#print("-----------------")
	prob = {}
	# calculate prob.
	for path in cells_paths:
		if level >= len(path):
			c = '' # top
		else:
			c = path[level]
		
		if c not in prob:
			prob[c] = 1
		else:
			prob[c] += 1

	# calculate entropy
	h = 0
	for cluster in prob:
		#print("$$$$$ " + str(cluster) + "           " + str(level) + "     " + str(cells_paths))
		new_cluster = [path for path in cells_paths if level < len(path) and path[level] == cluster];
		h += level_coeff*cluster_entropy(new_cluster, level+1)
		p = prob[cluster] / len(cells_paths);
		h -= p*math.log2(p)
	print(h)
	return h 


hdl_hier_list_s = [ cell.split('/')[0:-1] for cell in hdl_hier_list] # splitting 'a/b/c' to ['a', 'b', 'c']

hdl_hier_list_s = [
	['a','a','b1'],
	['a','a','b2'],
	['a','a','b2'],
	['b','b','b1'],
	['b','b','b3'],
	['b','b','b4']
]
cluster = [0,1,2,3,4,5]
cells_paths = [ hdl_hier_list_s[i] for i in cluster];
print("ENTROPY: " + str( cluster_entropy(cells_paths) ) )

assert False, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"






def detect_cluster(clusters_list, cell_path):
	for cluster in clusters_list:
		if re.match("^"+cluster+".*", '/' + cell_path):
			return cluster
	return "_NA_"


def dendrogram2tree(dendrogram, merge, expected_hier, gt_tree, tree, depth = 0):

	print(expected_hier)
	print('----------')

	cells = dendrogram.collapse_merge(merge)
	cells_hier = [ hdl_hier_list[cell] for cell in cells];
	cluster_entropy(cells_hier)

	assert False, "--------------------+++++++++++++==========================="


	split_points = dendrogram.split(len(expected_hier)*2, merge)



	histogram = {}
	for split in split_points:
		cells = dendrogram.collapse_merge(split)

		"""
		histogram[split] = {}#h:0 for h in expected_hier} #histogram[split]["_NA_"] = 0		
		for cell in cells:
			hdl_path = detect_cluster(expected_hier, hdl_hier_list[cell])
			if hdl_path not in histogram[split]:
				histogram[split][hdl_path] = 0

			histogram[split][hdl_path] += 1
		print("SPLIT - " + str(split) + ":")
		print(histogram[split])
		"""

		histo= {}
		# Option I - building expected tree according to gt_tree
		#for x in histogram[split].keys():
		#	print(x)
		#	expected_hier = gt_tree[x][0]
		#	print(expected_hier)

		
		# Option II - building expected tree according to split result
		# create histogram of all cells in split and filter out "tiny" clusters (% and/or abs size) THRESHOLD_p/THRESHOLD_abs
		for cell in cells:
			hdl_path = hdl_hier_list[cell].split('/')[0:-1]
			hdl_path = '/' + '/'.join(hdl_path)
			if hdl_path not in histo:
				histo[hdl_path] = 0
			histo[hdl_path] += 1

		THRESHOLD_p   = 0.02
		THRESHOLD_abs = 250			
		histo = { key:val for key,val in histo.items() if val/len(cells) > THRESHOLD_p or val >= THRESHOLD_abs}
		print(histo)
		for k in histo.keys():
			print(k)
		print("***************************************************")





my_tree = {}
#t.print()

expected_hier = []
cluster = ""
for hier in t.tree[cluster][0].keys():
	expected_hier.append(hier)


#if len(t.tree[cluster][1]):
#	if cluster == "":
#		cluster = '/'
#	expected_hier.append(cluster)

dendrogram2tree(dendo, dendo.root, expected_hier, t.tree, my_tree)
print(my_tree)

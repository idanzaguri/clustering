from igraph import *
import math
import re
"""
def flat_dendrogram(dendrogram, offset, cells):
	num_cells = dendrogram[-1][2]

	for i in range(2):
		if dendrogram[offset][i] < num_cells:
			cells.append(dendrogram[offset][i])
		else:
			flat_dendrogram(dendrogram, dendrogram[offset][i]-num_cells, cells)

	

def split_dendrogram(dendrogram, offset, tree):
	num_cells = dendrogram[-1][2]
	index = offset-num_cells;
	if index < 0:
		return
	num_of_sons = dendrogram[index][2]


	#print("split_dendrogram(%d) - %d"%(offset, num_of_sons))
	
	if num_of_sons <= 1000:
		l = [];
		flat_dendrogram(dendrogram, index , l)
		#print("%s" % (str(l)))
		tree[offset][1] = l;
		return


	exp_num_hier = round(math.log2(num_of_sons))
	#print(exp_num_hier)	
	assert exp_num_hier > 1 and exp_num_hier < 20, "Unexpected exp_num_hier"


	c = [dendrogram[index][0], dendrogram[index][1]]
	for i in range(exp_num_hier-2):
		m = max(c)
		c.remove(m)
		c.append(dendrogram[m-num_cells][0])
		c.append(dendrogram[m-num_cells][1])

	#print(c)
	
	tree[offset] = [[],[]]
	for i in c:
		if(i-num_cells < 0):
			tree[offset][1].append(i)			
		else:
			tree[offset][0].append(i)
			tree[i] = [[],[]]
			split_dendrogram(dendrogram, i, tree)
	
"""

def simplify_scikit_network_dendrogram(dendrogram):
	simplified_dendrogram = []
	for l in dendrogram:
		simplified_dendrogram.append([int(l[0]), int(l[1]), int(l[3])])
	return simplified_dendrogram


def network_dendrogram_2_igraph_dendrogram(dendrogram):
	simplified_dendrogram = []
	for l in dendrogram:
		simplified_dendrogram.append((int(l[0]), int(l[1])))
	return simplified_dendrogram





# CPP
# map <string, tuple<map <string, uint32_t> , vector<uint32_t>> >
# PYTHON
# { 'HIER' : [ { map of sons. son_hier : num_cells } , [list of cells] ]  }  
#

class HDLTree:
	tree = {}
	#SONS = 0
	#GATES = 1

	def __init__(self):
		self.tree[''] = [{},[]]
	
	
	def add_cell(self, cell, vid):
		m = re.split("/",cell)
		m.insert(0,'')	
		track = [ '/'.join(m[0:i+1]) for i in range(len(m)) ]

		for i in range(1,len(m)):
			parent = track[i-1]
			if i == len(m)-1: # cell
				self.tree[parent][1].append(vid)#track[i])
			else:	
				if track[i] not in self.tree:
					self.tree[parent][0][track[i]] = 0
					self.tree[track[i]] = [{},[]]
				self.tree[parent][0][track[i]] += 1


	def print(self, node = ''):
		depth = len(re.split("/",node))-1
		space = '    '*depth

		print(space+node+": ("+str(len(self.tree[node][1]))+")")
		for son in self.tree[node][0]:
			self.print(son)

		if node == '':
			print(len(self.tree.keys()))


	def num_of_cells(self, node):
		c = 0
		for son in self.tree[node][0]:
			c += self.tree[node][0][son]
		return c + len(self.tree[node][1])



class MyDendrogram:
	dendrogram = []

	num_dummy_vid = 0#len(hdl_hier_list)
	num_elements  = 0#len(hdl_hier_list)
	root          = 0#len(cpp_dendrogram)*2

	def __init__(self, dendrogram, num_dummy_vid = 0):
		if isinstance(dendrogram, VertexDendrogram):
			self.dendrogram = dendrogram.merges
		else:
			self.dendrogram = network_dendrogram_2_igraph_dendrogram(dendrogram)

		self.root = len(self.dendrogram)*2
		self.num_elements  = len(self.dendrogram)+1
		self.num_dummy_vid = num_dummy_vid
		

	def is_element(self, m):
		return (m < self.num_elements)

	def is_dummy_vertex(self, m):
		return (m < self.num_elements and m >= self.num_elements-self.num_dummy_vid)

	def get_merge_index(self, m):
		assert m >= self.num_elements, "ERROR: not a merge point!"
		return m - self.num_elements

	def split(self, n, merge):
		split_points = [merge]

		for i in range(n-1):
			m = split_points.pop( split_points.index(max(split_points)) )
			if self.is_element(m):
				print("ERROR")
				break # can't split anymore

			sub_merges = self.dendrogram[self.get_merge_index(m)]
			split_points.append(sub_merges[0])
			split_points.append(sub_merges[1])
		
		return split_points



	def collapse_merge(self, merge):
		assert not self.is_element(merge), "Can't collapse element!"
		cells = []
		merge_index = self.get_merge_index(merge)
		sub_merge = [0]*(merge_index+1)
		sub_merge[merge_index] = 1
	
		for index in range(merge_index,-1,-1):
			if sub_merge[index] == 0:
				continue
			for m in self.dendrogram[index]:
				if self.is_element(m):
					if not self.is_dummy_vertex(m):
						cells.append(m)
				else:
					sub_merge[self.get_merge_index(m)] = 1
		return cells
		#### recursive version - elegant but not applicable for python :\
		###if is_dummy_vertex(merge):
		###	print("DUMMY")
		###elif is_element(merge):
		###	print("ELEMNT")	
		###	cells.append(merge)
		###else:
		###	sub_merges = dendrogram[get_merge_index(merge)]
		###	print(sub_merges)		
		###	collapse_merge(dendrogram, sub_merges[0], cells)
		###	collapse_merge(dendrogram, sub_merges[1], cells)


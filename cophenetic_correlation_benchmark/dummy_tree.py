import sys
import random
import dendrogram_utils
from clustering_module import *

def num_sub_mod_prob(depth):
    opt = [0]*( int(depth/1.5) ** 3 )
    if depth <= 1:
        opt.extend(range(5,11))
    elif depth <= 3:
        opt.extend(range(1,8))
    elif depth <= 10: # Maximum tree depth
        opt.extend(range(1,6))
    return random.choice(opt)

def num_cells_prob(depth):
    if depth < 2:
        return random.randrange(0,100)
    if depth < 8:
        return random.randrange(500,1300)
    if depth < 15:
        return random.randrange(100, 250)
    return random.randrange(10,100)


def generate_hdl_hier_list(path, hdl_hier_list):

    depth = len(path.split('/'))
    
    for i in range(num_sub_mod_prob(depth)):
        generate_hdl_hier_list( path + "m" + str(i) + "/", hdl_hier_list)
    
    for i in range(num_cells_prob(depth)):
        hdl_hier_list.append(path + "c" + str(i))


random.seed(1628539009)

hdl_hier_list = []
generate_hdl_hier_list('', hdl_hier_list)
print(len(hdl_hier_list))
t = dendrogram_utils.HDLTree()

for vid,cell in enumerate(hdl_hier_list):
    t.add_cell(cell, vid)
    

mod_tree = t.tree
 
def print_tree(tree):
    for hier in sorted(tree.items()):
        print(hier[0])
#print_tree(mod_tree)


trimmed_hdl_hier_list = [h[1:] for h in hdl_hier_list]

swap_mode = 0
if len(sys.argv) >= 2:
	swap_mode = int(sys.argv[1])
zaguri_debug(mod_tree, hdl_hier_list, 4, 1, swap_mode)




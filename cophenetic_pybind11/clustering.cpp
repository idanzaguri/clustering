#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include <igraph.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <numeric>

#include <omp.h>

using namespace std;
namespace py = pybind11;

uint32_t MIN_CLUSTER_SIZE = 100;
uint32_t MAX_CLUSTER_SIZE = 1000;
uint32_t MAX_DEPTH        = 10;
uint32_t EXP_NUM_HIER = 0;

void set_MIN_CLUSTER_SIZE(uint32_t n) { MIN_CLUSTER_SIZE = n; }
void set_MAX_CLUSTER_SIZE(uint32_t n) { MAX_CLUSTER_SIZE = n; }
void set_MAX_DEPTH(uint32_t n) { MAX_DEPTH = n; }
void set_EXP_NUM_HIER(uint32_t n) { EXP_NUM_HIER = n; }


typedef map<string, tuple<map <string, uint32_t> , vector<uint32_t>> > hdl_tree_t;


pair<uint32_t, uint32_t> sorted_pair(uint32_t a, uint32_t b)
{
	return pair<uint32_t, uint32_t>( (a>=b)?a:b, (a>=b)?b:a );
}


/*
	         +---------------------------------------+
	root_id  | ( sons: next_id1, next_id2, cells:  ) |
	         +---------------------------------------|
	         |                                       |
	        ~~~                                     ~~~
	         |                                       |
	         +---------------------------------------+
	next_id2 | ( sons: , cells: xxx )                |
	         +---------------------------------------|
	         |                                       |
	        ~~~                                     ~~~
	         |                                       |
	         +---------------------------------------+
	next_id1 | ( sons: , cells: xxx )                |
	         +---------------------------------------|
	
*/

class HDLTree {
	public:
		// python compatible data structure		
		typedef map <uint32_t, tuple< vector<uint32_t>, vector<uint32_t>>> PC_type;

	private:
		uint32_t root;
		uint32_t num_elements;
		PC_type tree;
		
		map <uint32_t, uint32_t> cell2node;
		map <pair<uint32_t, uint32_t>, int> n2n_dist;

		map <uint32_t, vector<uint32_t>> nodes_track_map;
		void calc_nodes_track(uint32_t node, vector <uint32_t> track)
		{
			track.push_back(node);
		//	cout << "DEBUG!!!!" << node << " " << std::get<0>(tree.at(node)).size() << " " << track.size() << endl;
			nodes_track_map.at(node) = track;
			vector<uint32_t> sons = std::get<0>(tree.at(node));
			for(size_t i=0; i<sons.size(); i++)
			{
		//		cout << i << "/" << sons.size() << endl;
				calc_nodes_track(sons[i], track);
			}
		}

		void calc_nodes_track()
		{
			for(PC_type::iterator it = tree.begin(); it != tree.end(); ++it)
				nodes_track_map[it->first];
			calc_nodes_track(root, vector <uint32_t>());

			for(PC_type::iterator it1 = tree.begin(); it1 != tree.end(); ++it1) {
				for(PC_type::iterator it2 = tree.begin(); it2 != tree.end(); ++it2) {
					if(it1->first < it2->first) continue;

					vector<uint32_t> t1 = nodes_track_map.at(it1->first);
					vector<uint32_t> t2 = nodes_track_map.at(it2->first);

					size_t lca;
					assert(t1.size() && t2.size());
					for(lca=0; lca<t1.size() && lca<t2.size() && t1[lca]==t2[lca]; ++lca) {}
					n2n_dist[ sorted_pair(it1->first, it2->first) ] = t1.size() + t2.size() - 2*(lca);
				}
			}
		}


		int node_dist(uint32_t node1, uint32_t node2) { 
			return n2n_dist.at(sorted_pair(node1, node2)); 
		}
	public:

		void add_node(uint32_t node)
		{
			PC_type::iterator it = tree.find(node);
			assert (it == tree.end()); // Duplicate addition
			tree[node];
		}


		auto& get_cells(uint32_t node) { return std::get<1>(tree.at(node));}

		void add_son (uint32_t node, uint32_t x) { std::get<0>(tree.at(node)).push_back(x); }
		void add_cell(uint32_t node, uint32_t x) { std::get<1>(tree.at(node)).push_back(x); }

		void sort()
		{
			for(PC_type::iterator it = tree.begin(); it != tree.end(); ++it)
			{
				vector<uint32_t> *sons = &(std::get<0>(it->second));
				vector<uint32_t> *cells = &(std::get<1>(it->second));
				std::sort (sons->begin(), sons->end());
				std::sort (cells->begin(), cells->end());
			}
		}

		void set_root(uint32_t root) { this->root = root;}
		void set_num_elements(uint32_t num_elements) { this->num_elements = num_elements;}

		vector<uint32_t> get_distance_vector(uint32_t cell_id);

		void initilaize_speedup_maps();

		uint32_t sum=0;

		void print(uint32_t node, int depth)
		{
			for(int i=0; i<depth; i++) cout << "\t";
			
			//assert (std::get<1>(tree.at(node)).size() >= MIN_CLUSTER_SIZE);

			cout << "(" << node << ")    " << std::get<1>(tree.at(node)).size() << endl;
			sum+=std::get<1>(tree.at(node)).size();
		
			for(size_t i=0; i<(std::get<0>(tree.at(node))).size(); ++i)
			{
				print((std::get<0>(tree.at(node)))[i], depth+1);
			}
		}

		void print() {sum =0; print(root,0); cout << "NUM_OF_CELLS:" << sum << endl; cout << "NUM_OF_NODES: " <<  tree.size() << endl;}
		uint32_t size() { return num_elements; }


};







void HDLTree::initilaize_speedup_maps()
{
	for(PC_type::const_iterator it = tree.begin(); it != tree.end(); ++it) {
		vector<uint32_t> cells = std::get<1>(it->second);
		for(size_t i=0; i<cells.size(); ++i)
			 cell2node[ cells[i] ] = it->first;
	}
	calc_nodes_track();
}






vector<uint32_t> HDLTree::get_distance_vector(uint32_t cell_id)
{
	// initialize distance vector
	vector<uint32_t> dist_vec(num_elements-(cell_id+1), 1e6);
	uint32_t node = cell2node.at(cell_id);

	for(uint32_t i=0; i<num_elements-(cell_id+1); i++) {
		dist_vec[i] = node_dist(node, cell2node.at(cell_id+1+i)); 
	}
	
	return dist_vec;
}












class Dendrogram {
	public:
		// python compatible data structure
		typedef vector< tuple<uint32_t, uint32_t, uint32_t> > PC_type;
		typedef vector< tuple<uint32_t, uint32_t> > PC_igraph_type;

	private:
		uint32_t num_elements;
		uint32_t num_dummy_vids;

		uint32_t root;

		PC_type dendrogram;
	
		uint32_t add_num_of_leafs(size_t merge);

		inline uint32_t get_merge_index(size_t merge) { return merge-num_elements; }
		inline bool is_element(size_t merge) { return (merge < num_elements); }

		inline bool  is_dummy_vertex(size_t merge) { return (merge >= num_elements-num_dummy_vids && merge < num_elements); }

		uint32_t get_0(size_t merge)        { return std::get<0>(dendrogram[ get_merge_index(merge) ]); }
		uint32_t get_1(size_t merge)        { return std::get<1>(dendrogram[ get_merge_index(merge) ]); }
		uint32_t get_num_sons(size_t merge) { return std::get<2>(dendrogram[ get_merge_index(merge) ]); }


		void convert_to_tree(size_t merge, HDLTree &tree, uint32_t depth);
	public:
		Dendrogram(PC_type dendrogram) { assert(false); this->dendrogram = dendrogram; num_elements = dendrogram.size()+1; root = dendrogram.size()*2;}


		Dendrogram(PC_igraph_type dendrogram, uint32_t _num_dummy_vids = 0) {
			num_dummy_vids = _num_dummy_vids;

			for(size_t i = 0; i<dendrogram.size(); i++) {
				this->dendrogram.push_back( std::make_tuple( get<0>(dendrogram[i]), get<1>(dendrogram[i]), -1 ) );
			}
			num_elements = dendrogram.size()+1;
			
			root = dendrogram.size()*2;
			assert (add_num_of_leafs(root) == num_elements);
		}

		void collapse_merge(size_t merge, vector<uint32_t> &elements);
		HDLTree convert_to_tree();


		HDLTree convert_to_tree_NEW(hdl_tree_t gt_tree);		
		void convert_to_tree_NEW(size_t merge, string cluster, hdl_tree_t gt_tree, HDLTree &tree);

};


uint32_t Dendrogram::add_num_of_leafs(size_t merge)
{
	if(is_element(merge)) return 1;
	uint32_t n = 
		add_num_of_leafs( get_0(merge) ) +
		add_num_of_leafs( get_1(merge) );

	std::get<2>(dendrogram[ get_merge_index(merge) ]) = n;
	return n;
}






void Dendrogram::collapse_merge(size_t merge, vector<uint32_t> &elements)
{
	if (is_element(get_0(merge))) elements.push_back(get_0(merge));
	else collapse_merge(get_0(merge), elements);

	if (is_element(get_1(merge))) elements.push_back(get_1(merge));
	else collapse_merge(get_1(merge), elements);
}


void Dendrogram::convert_to_tree(size_t merge, HDLTree &tree, uint32_t depth)
{
	if (is_element(merge)) {
		assert(false); // ERROR
		return;
	}

	tree.add_node(merge);	

	//assert (get_num_sons(merge) > MIN_CLUSTER_SIZE); FIXME: 
	
	if(get_num_sons(merge) <= MAX_CLUSTER_SIZE || depth >= MAX_DEPTH) {
		collapse_merge(merge, tree.get_cells(merge));
		return;
	}

	uint32_t exp_num_hier = (EXP_NUM_HIER==0)? round(log2(get_num_sons(merge))) : EXP_NUM_HIER;
	assert(exp_num_hier > 1 && exp_num_hier < 50); // Unexpected exp_num_hier




	vector <uint32_t> c;
	c.push_back( get_0(merge) ); std::push_heap(c.begin(), c.end());
	c.push_back( get_1(merge) ); std::push_heap(c.begin(), c.end());

	for(size_t i=0; i < exp_num_hier-2; i++) {

		std::pop_heap(c.begin(), c.end()); // moves the largest to the end
		uint32_t m = c.back();
		c.pop_back();  // actually removes the largest element

		//std::cout << "C: "; for(size_t  jj=0; jj <c.size(); jj++) std::cout << c[jj] << ", "; std::cout << std::endl;

		if(is_element(m)) {
			tree.add_cell(merge, m);
			i--;
		} else if(get_num_sons(m) < MIN_CLUSTER_SIZE) {	
			collapse_merge(m, tree.get_cells(merge));
			i--;
		} else {
			c.push_back( get_0(m) ); std::push_heap(c.begin(), c.end());
			c.push_back( get_1(m) ); std::push_heap(c.begin(), c.end());
		}
	}
	

	for(size_t i=0; i < c.size(); i++) {
		if(is_element(c[i])) {
			tree.add_cell(merge, c[i]);
		}
		else {
			if(get_num_sons(c[i]) < MIN_CLUSTER_SIZE) {
				collapse_merge(c[i], tree.get_cells(merge));
			} else {
				tree.add_son(merge, c[i]);
				convert_to_tree(c[i], tree,depth+1);
			}
		}	
	}

}


HDLTree Dendrogram::convert_to_tree()
{
	HDLTree tree;
	tree.set_root(root);
	tree.set_num_elements(num_elements);
	convert_to_tree(root, tree, 0);
	tree.sort();
	return tree;
}



int hdl_hierarchy_distance(const char *cell1, const char *cell2)
{
	int s1 = strlen(cell1);
	int s2 = strlen(cell2);
	int size = (s1<s2)? s1 : s2;
	int hier_counter=0;
	for(int i=0; i<size; i++)
	{
		if(cell1[i] == cell2[i]) continue;
		for(int j=i; j<s1; j++) if(cell1[j] == '/') ++hier_counter;
		for(int j=i; j<s2; j++) if(cell2[j] == '/') ++hier_counter;
		break;
	}
	return hier_counter;
}








double cophenetic_correlation(Dendrogram dendrogram, vector<string> hdl_hierarchies, int request_nthreads, int debug)
{
	
	HDLTree tree = dendrogram.convert_to_tree();
	tree.print();
	tree.initilaize_speedup_maps();
	int nthreads;

	double hdl_hierarchy_distance_mean = 0;
	uint64_t hdl_vec_size = 0;

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_hdl_hierarchy_distance_mean = 0;
		uint64_t l_hdl_vec_size = 0;

		for (long unsigned int i = id; i < hdl_hierarchies.size(); i+=nthrds) {
			for (long unsigned int j = i+1; j < hdl_hierarchies.size(); j++) {
				l_hdl_hierarchy_distance_mean += hdl_hierarchy_distance(hdl_hierarchies[i].c_str(), hdl_hierarchies[j].c_str());
				++l_hdl_vec_size;
			}
			
			if(debug && ((i & 0xFFF) == 0))
				cout << id << ": " << i << "/" << hdl_hierarchies.size() << endl;
		}
		
		#pragma omp critical
		{
			hdl_hierarchy_distance_mean += l_hdl_hierarchy_distance_mean;
			hdl_vec_size += l_hdl_vec_size;
		}

	}
	hdl_hierarchy_distance_mean /= hdl_vec_size;

	if(debug) {
		cout << "#threads: " << nthreads <<  endl;
		cout << "hdl_hierarchies.size(): " << hdl_hierarchies.size() <<  endl;	
		cout << "hdl_vec_size: " << hdl_vec_size << endl;	
		cout << "hdl_hierarchy_distance_mean: " << hdl_hierarchy_distance_mean << endl;
		cout << "-------------------------------------------------" << endl;
	}


	double tree_distances_mean = 0;
	uint64_t tree_distances_vec_size = 0;

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_tree_distances_mean = 0;
		uint64_t l_tree_distances_vec_size = 0;

		for (long unsigned int i = id; i < tree.size(); i+=nthrds) {
			vector<uint32_t> d = tree.get_distance_vector(i);
			l_tree_distances_mean += accumulate(d.begin(), d.end(), 0);
			l_tree_distances_vec_size += d.size();
			if(debug && ((i & 0xFFF)  == 0))
				cout << id << ": " << i << "/" << tree.size() << endl;
		}

		#pragma omp critical
		{
			tree_distances_mean     += l_tree_distances_mean;
			tree_distances_vec_size += l_tree_distances_vec_size;

		}

	}
	tree_distances_mean /= tree_distances_vec_size;

	if(debug) {
		cout << "#threads: " << nthreads << endl;
		cout << "tree.size(): " << tree.size() << endl;
		cout << "tree_distances_vec_size: " << tree_distances_vec_size << endl;
		cout << "tree_distances_mean: " << tree_distances_mean << endl;
		cout << "-------------------------------------------------" << endl;
	}


	double c_numerator, c_denominator1, c_denominator2;

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{

		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_c_numerator = 0, l_c_denominator1 = 0, l_c_denominator2 = 0;

		for(long unsigned int i=id; i<tree.size(); i+=nthreads) {
			vector<uint32_t> tree_vec = tree.get_distance_vector(i);
			vector<uint32_t> hdl_vec;
			for (long unsigned int j = i+1; j < hdl_hierarchies.size(); j++)
				hdl_vec.push_back(hdl_hierarchy_distance(hdl_hierarchies[i].c_str(), hdl_hierarchies[j].c_str()));

			assert (tree_vec.size() == hdl_vec.size());

			for(long unsigned int j=0; j<tree_vec.size(); j++) {
				double hdl_tmp = hdl_vec[j]-hdl_hierarchy_distance_mean;
				double den_tmp = tree_vec[j]-tree_distances_mean;
				l_c_numerator += hdl_tmp*den_tmp;
				l_c_denominator1 += pow(hdl_tmp,2);
				l_c_denominator2 += pow(den_tmp,2);
			}
	
			if(debug && ((i & 0xFFF) == 0))
				cout << id << ": " << i << "/" << tree.size() << endl;
		}

		#pragma omp critical
		{
			c_numerator += l_c_numerator;
			c_denominator1 += l_c_denominator1;
			c_denominator2 += l_c_denominator2;
		}
	}
	double cophenetic_factor = c_numerator / sqrt(c_denominator1*c_denominator2);

	if(debug) {
		cout << "#threads: " <<  nthreads <<  endl;
		cout << "RESULT: " << cophenetic_factor << endl;
	}

	return cophenetic_factor;
}





double cophenetic_correlation(Dendrogram::PC_type d, vector<string> hdl_hierarchies, int nthreads, int debug)
{
	return cophenetic_correlation(Dendrogram(d), hdl_hierarchies, nthreads, debug);
}

double cophenetic_correlation(Dendrogram::PC_igraph_type d, vector<string> hdl_hierarchies, int nthreads, int debug)
{
	return cophenetic_correlation(Dendrogram(d), hdl_hierarchies, nthreads, debug);
}









































































void Dendrogram::convert_to_tree_NEW(size_t merge, string cluster, hdl_tree_t gt_tree, HDLTree &tree)
{



	if (is_element(merge)) {
		assert(false); // ERROR is_element
		return;
	}
	if (is_dummy_vertex(merge)) {
		assert(false); // ERROR is_dummy_vertex
		return;
	}



	map <string, uint32_t> sons = std::get<0>(gt_tree.at(cluster));
	for(map <string, uint32_t> ::iterator it = sons.begin(); it != sons.end(); ++it)
		cout << it->first << endl;



	tree.add_node(merge);	


	uint32_t exp_num_hier = sons.size();
	assert(exp_num_hier > 1 && exp_num_hier < 50); // Unexpected exp_num_hier


	vector <uint32_t> c;
	c.push_back( get_0(merge) ); std::push_heap(c.begin(), c.end());
	c.push_back( get_1(merge) ); std::push_heap(c.begin(), c.end());

	for(size_t i=0; i < exp_num_hier-2; i++) {
		std::pop_heap(c.begin(), c.end()); // moves the largest to the end
		uint32_t m = c.back();
		c.pop_back();  // actually removes the largest element

		//std::cout << "C: "; for(size_t  jj=0; jj <c.size(); jj++) std::cout << c[jj] << ", "; std::cout << std::endl;

		if(is_dummy_vertex(m))
		{
			cout << "ERROR: found dummy vertex!!!\n";
			assert(false);
		} 
		else if(is_element(m)) 
		{
			cout << "ERROR: found element!!!\n";
			assert(false);
			//tree.add_cell(merge, m);
			//i--;
		}
		else
		{
			c.push_back( get_0(m) ); std::push_heap(c.begin(), c.end());
			c.push_back( get_1(m) ); std::push_heap(c.begin(), c.end());
		}
	}

	
/*
	for(size_t i=0; i < c.size(); i++) {
		if(is_element(c[i])) {
			tree.add_cell(merge, c[i]);
		}
		else {
			if(get_num_sons(c[i]) < MIN_CLUSTER_SIZE) {
				collapse_merge(c[i], tree.get_cells(merge));
			} else {
				tree.add_son(merge, c[i]);
				convert_to_tree(c[i], tree,depth+1);
			}
		}	
	}

*/


}



HDLTree Dendrogram::convert_to_tree_NEW(hdl_tree_t gt_tree)
{
	HDLTree tree;
	tree.set_root(root);
	tree.set_num_elements(num_elements);
	convert_to_tree_NEW(root, "", gt_tree, tree);
	tree.sort();
	return tree;	
}	



//typedef map<string, tuple<map <string, uint32_t> , vector<uint32_t>> > hdl_tree_t;
void print_hdl_tree(Dendrogram::PC_igraph_type d, hdl_tree_t gt_tree, uint32_t num_dummy_vids, vector<string> hdl_hierarchies)
{
	Dendrogram dendrogram(d, num_dummy_vids);
	dendrogram.convert_to_tree_NEW(gt_tree);
//	for(hdl_tree_t::iterator it = gt_tree.begin(); it != gt_tree.end(); ++it)
//		cout << it->first << endl;
}







PYBIND11_MODULE(clustering_module, handle) {
	handle.doc() = "clustering_module... Idan Zaguri";
	handle.def("cophenetic_correlation", py::overload_cast<Dendrogram::PC_type, vector<string>, int, int>(&cophenetic_correlation));
	handle.def("cophenetic_correlation", py::overload_cast<Dendrogram::PC_igraph_type, vector<string>, int, int>(&cophenetic_correlation));

	handle.def("set_MIN_CLUSTER_SIZE", &set_MIN_CLUSTER_SIZE);
	handle.def("set_MAX_CLUSTER_SIZE", &set_MAX_CLUSTER_SIZE);
	handle.def("set_MAX_DEPTH", &set_MAX_DEPTH);
	handle.def("set_EXP_NUM_HIER", &set_EXP_NUM_HIER);

	handle.def("print_hdl_tree", &print_hdl_tree);
}

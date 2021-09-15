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
#include <random>

#include <omp.h>

using namespace std;
namespace py = pybind11;




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
		typedef map <uint32_t, tuple< vector<uint32_t>       , vector<uint32_t>> > PC_type;
		typedef map <string  , tuple< map <string, uint32_t> , vector<uint32_t>> > PC_type_new;
		map < string, uint32_t > hier2id;

		HDLTree(PC_type_new py_tree)
		{
			num_elements = 0;

			uint32_t unique_id = 0;

			// for each node in tree
			for (PC_type_new::iterator it = py_tree.begin(); it != py_tree.end(); ++it)
			{
				string hier = it->first;

				if(hier2id.find(hier) == hier2id.end()) {
					hier2id.insert(pair<string, uint32_t>(hier, unique_id));
					tree[unique_id++];
				}
				uint32_t node = hier2id[hier];

				map <string, uint32_t> children = std::get<0>(it->second);
				vector <uint32_t>      cells    = std::get<1>(it->second);


				// for each child in node convert child path to unique_id
				for (map <string, uint32_t>::iterator it1 = children.begin(); it1 != children.end(); ++it1)
				{
					string child = it1->first;
					if(hier2id.find(child) == hier2id.end()) {
						hier2id.insert(pair<string, uint32_t>(child, unique_id));
						tree[unique_id++];
					}
					std::get<0>(tree.at(node)).push_back( hier2id[child] );
				}

				std::get<1>(tree.at(node)) = cells; // copy cells list from py_tree to tree
				num_elements += cells.size();
			}

			root = hier2id["/"];
		}



		void remove_node_from_parent(uint32_t node)
		{
			// for each node in tree find if it node's parent
			for (PC_type::iterator it = tree.begin(); it != tree.end(); ++it)
			{
				uint32_t parent = it->first;
				cout << parent << endl;
				vector<uint32_t> &children = std::get<0>(tree.at(parent));

				uint32_t pre_size = children.size();
				children.erase( std::remove(children.begin(), children.end(), node) , children.end());
				if(pre_size != children.size()) return;
			}
			assert(false); // can't locate parent!
		}

		void move_hier(uint32_t src, uint32_t dst)
		{
			remove_node_from_parent(src);
			std::get<0>(tree.at(dst)).push_back(src);
		}


		void flat_hier(uint32_t node, vector<uint32_t> &flat_cells, int keep = 1)
		{
			vector<uint32_t> children = std::get<0>(tree.at(node));
			vector<uint32_t> cells    = std::get<1>(tree.at(node));
			flat_cells.insert(std::end(flat_cells), std::begin(cells), std::end(cells));
				
			for (auto i = children.begin(); i != children.end(); ++i) {
				flat_hier(*i, flat_cells, 0);
			}
			if(!keep)
				tree.erase(node);

		}

		void flat_hier(uint32_t node)
		{
			vector<uint32_t> flat_cells;
			flat_hier(node, flat_cells);
			std::get<0>(tree.at(node)).clear();
			std::get<1>(tree.at(node)) = flat_cells;
		}
		


		void keep_firts_level_only()
		{
			vector<uint32_t> level0 = std::get<0>(tree.at(0));
			for(auto i = level0.begin(); i != level0.end(); ++i) {
				flat_hier(*i);
			}
		}




		void shuffle_clean_tree(uint32_t node, vector<uint32_t> &cells)
		{
			vector<uint32_t>  children   = std::get<0>(tree.at(node));
			vector<uint32_t> &node_cells = std::get<1>(tree.at(node));
			cells.insert(std::end(cells), std::begin(node_cells), std::end(node_cells));

			uint32_t node_size = node_cells.size();
			node_cells.clear();
			node_cells.push_back(node_size);

			for (auto i = children.begin(); i != children.end(); ++i) {
				shuffle_clean_tree(*i, cells);
			}
		}


		void shuffle_repopulate(uint32_t node, vector<uint32_t> &cells)
		{
			vector<uint32_t>  children   = std::get<0>(tree.at(node));
			vector<uint32_t> &node_cells = std::get<1>(tree.at(node));

			uint32_t node_size = node_cells.back();
			node_cells.pop_back();
			assert(node_cells.size() == 0);

			auto it = std::next(cells.begin(), node_size);	
			std::move(cells.begin(), it, std::back_inserter(node_cells));
			cells.erase(cells.begin(), it );

			for (auto i = children.begin(); i != children.end(); ++i) {
				shuffle_repopulate(*i, cells);
			}
		}


		void shuffle()
		{
			vector<uint32_t> cells;
			shuffle_clean_tree(0, cells); // get all cells in tree
			auto rng = std::default_random_engine {};
			std::shuffle(std::begin(cells), std::end(cells), rng);
			shuffle_repopulate(0, cells);
			assert (cells.size() == 0);
		}

		void spread_hier(uint32_t node)
		{
			vector<uint32_t> children = std::get<0>(tree.at(node));
			vector<uint32_t> cells    = std::get<1>(tree.at(node));
			assert(children.size() == 0); // need to run flat_hier before
			std::get<1>(tree.at(node)).clear();

			assert( (std::get<0>(tree.at(node))).size() == 0 && (std::get<1>(tree.at(node))).size() == 0 ); // should be empty node				
			tree.erase(node);

			while(cells.size()) {
				PC_type::iterator it;
				for( it = tree.begin(); it != tree.end() && cells.size() ; it++) {
					(std::get<1>(it->second)).push_back(cells.back());
					cells.pop_back();
				}
			}
			remove_node_from_parent(node);
		}



	private:
		uint32_t root;
		uint32_t num_elements;
		PC_type tree;

		map <uint32_t, uint32_t> cell2node;
		map <pair<uint32_t, uint32_t>, int> n2n_dist;

		map <uint32_t, vector<uint32_t>> nodes_track_map;
		void calc_nodes_track(uint32_t node, vector <uint32_t> track)
		{
			//cout << track.size() << endl;
			track.push_back(node);
			//cout << "DEBUG!!!!" << node << " " << std::get<0>(tree.at(node)).size() << " " << track.size() << endl;
			nodes_track_map.at(node) = track;
			vector<uint32_t> sons = std::get<0>(tree.at(node));
			for(size_t i=0; i<sons.size(); i++)
				calc_nodes_track(sons[i], track);
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


double cophenetic_correlation(HDLTree tree, vector<string> hdl_hierarchies, int request_nthreads, int debug, int swap)
{
	if(swap == 1) {
		tree.move_hier(8, 174);
	}
	if(swap == 2) {
		tree.flat_hier(6);
	}
	if(swap == 22) {
		tree.flat_hier(6);
		tree.spread_hier(6);
	}
	if(swap == 3) {
		tree.flat_hier(0);
	}
	if(swap == 4) {
		tree.shuffle();
	}
	if(swap == 5) {
		tree.keep_firts_level_only();
	}

	tree.sort();
	tree.initilaize_speedup_maps();

	if(debug) tree.print();

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
	cout << c_numerator << endl;
	cout << c_denominator1 << endl;
	cout << c_denominator2 << endl;

	double cophenetic_factor = c_numerator / sqrt(c_denominator1*c_denominator2);
		
	if(isnan(cophenetic_factor)) {
		cophenetic_factor = 0;
		cout << "NaN -> 0" << endl;		
	}
	if(debug) {
		cout << "#threads: " <<  nthreads <<  endl;
		cout << "RESULT: " << cophenetic_factor << endl;
	}

	return cophenetic_factor;
}


double zaguri_debug(HDLTree::PC_type_new py_tree, vector<string> hdl_hierarchies, int request_nthreads, int debug, int swap)
{
	return cophenetic_correlation(HDLTree(py_tree), hdl_hierarchies, request_nthreads, debug, swap);
}

PYBIND11_MODULE(clustering_module, handle) {
	handle.doc() = "clustering_module... Idan Zaguri";
	handle.def("zaguri_debug", &zaguri_debug);
}

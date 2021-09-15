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

#include <omp.h>

using namespace std;
namespace py = pybind11;



class Dendrogram {
	private:
		vector< vector<uint32_t> > dendrogram;
		uint32_t num_of_gates;

	public:
		Dendrogram(vector< vector<uint32_t> > _dendrogram) : dendrogram (_dendrogram) { num_of_gates = _dendrogram.size()+1;}

		inline uint32_t find_index_location(uint32_t index, uint32_t start);

		void calc_dist_vec(vector<uint32_t> &dist_vec, uint32_t src_index, uint32_t index, uint32_t loc, uint32_t depth, int direction_up);
		vector<uint32_t> get_distance_vector(uint32_t index);
};

inline uint32_t Dendrogram::find_index_location(uint32_t index, uint32_t start = 0)
{

	if(/*is not gate*/ index > num_of_gates) start = index-num_of_gates-1;

	for (long unsigned int i = start; i < dendrogram.size(); i++)
		if(dendrogram[i][0] == index || dendrogram[i][1] == index)
			return i;
	cout << "ERROR [find_index_location]: can't find requested index " << index <<"----"  << num_of_gates << endl;
	return 0;
}


void Dendrogram::calc_dist_vec(vector<uint32_t> &dist_vec, uint32_t src_index, uint32_t index, uint32_t loc, uint32_t depth, int direction_up = 1)
{
	//DEBUG cout << "index:" << index << "   location:" << loc << "    depth:" << depth << endl;

	//if(depth > 10) return;
	//if(index >= dendrogram.size()*2)
	//	return;

	for(int i=0; i<2; i++) {
		uint32_t neighbor = dendrogram[loc][i];
		if(neighbor == index || neighbor < src_index) continue;

		//DEBUG cout << neighbor;

		if(/*is gate*/ neighbor < num_of_gates) {
			//DEBUG cout << " -neighbor is gate" << endl;
			if(dist_vec[neighbor-src_index-1] != 1e6) cout << "ERROR" << endl;
			dist_vec[neighbor-src_index-1] = depth;
		} else {
			//DEBUG cout << " -neighbor is not gate" << endl;
			calc_dist_vec(dist_vec, src_index, neighbor, neighbor-dendrogram.size()-1, depth+1, 0);
		}
	}

	if(!direction_up) return;

	// go up to next merge
	index = loc+1+dendrogram.size();
	if(index >= dendrogram.size()*2)
		return;
	//DEBUG cout << index <<" going up!" << endl;
	loc = find_index_location(index);
	calc_dist_vec(dist_vec, src_index, index, loc, depth+1);
}


vector<uint32_t> Dendrogram::get_distance_vector(uint32_t index)
{
	// initialize distance vector
	vector<uint32_t> dist_vec(num_of_gates-index-1,1e6);

	calc_dist_vec(dist_vec, index, index ,find_index_location(index), 0);

	// check distance vector
	for(long unsigned int i = 0; i <= dist_vec.size(); i++) {
		if(dist_vec[i] == 1e6) {
			cout << "ERROR: (" << index << ") " << i << endl;
			break;
		}
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




double cophenetic_correlation(vector< vector<uint32_t> > _dendrogram, vector < string > hdl, int request_nthreads=1) {

	Dendrogram dendrogram(_dendrogram);

	int nthreads;

	double dendrogram_distance_mean = 0;
	uint64_t dendrogram_vec_size = 0;

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_dendrogram_distance_mean = 0;
		uint64_t l_dendrogram_vec_size = 0;

		//for(long unsigned int i=id; i<_dendrogram.size(); i+=nthrds) {
		for(long unsigned int i=id; i<_dendrogram.size(); i+=nthrds) {
			vector<uint32_t> vec = dendrogram.get_distance_vector(i);
			for (long unsigned int j = 0; j < vec.size(); j++) {
				l_dendrogram_distance_mean += vec[j];
				++l_dendrogram_vec_size;
			}
			if( (i & 0xFFF)  == 0)
				cout << id << ": " << i << "/" << _dendrogram.size() << endl;
		}
		#pragma omp critical
		{
			dendrogram_distance_mean += l_dendrogram_distance_mean;
			dendrogram_vec_size += l_dendrogram_vec_size;
		}

	}

	cout << "#threads: " <<  nthreads <<  endl;
	dendrogram_distance_mean /= dendrogram_vec_size;
	cout << "dendrogram_distance_mean: " <<  dendrogram_distance_mean <<  endl;


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

		for (long unsigned int i = id; i < hdl.size(); i+=nthrds) {
			for (long unsigned int j = i+1; j < hdl.size(); j++) {
				l_hdl_hierarchy_distance_mean += hdl_hierarchy_distance(hdl[i].c_str(), hdl[j].c_str());
				++l_hdl_vec_size;
			}
			if( (i & 0xFFF)  == 0)
				cout << id << ": " << i << "/" << hdl.size() << endl;
		}

		#pragma omp critical
		{
			hdl_hierarchy_distance_mean += l_hdl_hierarchy_distance_mean;
			hdl_vec_size += l_hdl_vec_size;

		}

	}

	cout << "#threads: " <<  nthreads <<  endl;
	hdl_hierarchy_distance_mean /= hdl_vec_size;
	cout << "hdl_hierarchy_distance_mean: " <<  hdl_hierarchy_distance_mean <<  endl;


	if(dendrogram_vec_size != hdl_vec_size) {
		cout << "ERROR: " << dendrogram_vec_size << " != " << hdl_vec_size << endl;
		return -1;
	}



	double c_numerator, c_denominator1, c_denominator2;

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{

		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_c_numerator = 0, l_c_denominator1 = 0, l_c_denominator2 = 0;

		for(long unsigned int i=id; i<_dendrogram.size(); i+=nthreads) {
			vector<uint32_t> dendrogram_vec = dendrogram.get_distance_vector(i);
			vector<uint32_t> hdl_vec;
			for (long unsigned int j = i+1; j < hdl.size(); j++)
				hdl_vec.push_back(hdl_hierarchy_distance(hdl[i].c_str(), hdl[j].c_str()));

			for(long unsigned int j=0; j<dendrogram_vec.size(); j++) {
				double hdl_tmp = hdl_vec[j]-hdl_hierarchy_distance_mean;
				double den_tmp = dendrogram_vec[j]-dendrogram_distance_mean;

				l_c_numerator += hdl_tmp*den_tmp;
				l_c_denominator1 += pow(hdl_tmp,2);
				l_c_denominator2 += pow(den_tmp,2);
			}
			if( (i & 0xFFF)  == 0)
				cout << id << ": " << i << "/" << _dendrogram.size() << endl;
		}

		#pragma omp critical
		{
			c_numerator += l_c_numerator;
			c_denominator1 += l_c_denominator1;
			c_denominator2 += l_c_denominator2;
		}
	}
	cout << "#threads: " <<  nthreads <<  endl;

	double cophenetic_factor = c_numerator / sqrt(c_denominator1*c_denominator2);
	cout << "RESULT: " << cophenetic_factor << endl;
	return cophenetic_factor;
}




#include <cassert>












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

typedef map <int, pair< vector<uint32_t>, vector<uint32_t>>> hdl_tree_t;




typedef vector < string > hdl_hier_t;
typedef vector< tuple<uint32_t, uint32_t, uint32_t> > dendrogram_t;
typedef vector< tuple<uint32_t, uint32_t> > igraph_dendrogram_merges_t;







void flat_dendrogram(const dendrogram_t &dendrogram, unsigned long int offset, vector<uint32_t> &cells)
{
	uint32_t num_leafs = dendrogram.size()+1;

	if (std::get<0>(dendrogram[offset]) < num_leafs)
		cells.push_back(std::get<0>(dendrogram[offset]));
	else
		flat_dendrogram(dendrogram, std::get<0>(dendrogram[offset])-num_leafs, cells);

	if (std::get<1>(dendrogram[offset]) < num_leafs)
		cells.push_back(std::get<1>(dendrogram[offset]));
	else
		flat_dendrogram(dendrogram, std::get<1>(dendrogram[offset])-num_leafs, cells);
}

void dendrogram_to_tree(const dendrogram_t &dendrogram, unsigned long int offset, hdl_tree_t &tree)
{
	uint32_t num_leafs = dendrogram.size()+1;
	//uint32_t num_cells = std::get<2>(dendrogram.back());
	long int index = offset-num_leafs;

	if (index < 0)
		return;

	uint32_t num_of_sons = std::get<2>(dendrogram[index]);

	if(num_of_sons <= 1000) {
		flat_dendrogram(dendrogram, index , std::get<1>(tree[offset]));
		return;
	}

	uint32_t exp_num_hier = round(log2(num_of_sons));

	assert(exp_num_hier > 1 && exp_num_hier < 20); // Unexpected exp_num_hier


	vector <uint32_t> c;
	c.push_back( std::get<0>(dendrogram[index]) ); std::push_heap(c.begin(), c.end());
	c.push_back( std::get<1>(dendrogram[index]) ); std::push_heap(c.begin(), c.end());

	for(unsigned long int i=0; i < exp_num_hier-2; i++) {

		std::pop_heap(c.begin(), c.end()); // moves the largest to the end
		uint32_t m = c.back();
		c.pop_back();  // actually removes the largest element

		c.push_back( std::get<0>(dendrogram[m-num_leafs]) ); std::push_heap(c.begin(), c.end());
		c.push_back( std::get<1>(dendrogram[m-num_leafs]) ); std::push_heap(c.begin(), c.end());
	}


	for(unsigned long int i=0; i < c.size(); i++) {
		if(c[i]-num_leafs < 0) {
			std::get<1>(tree[offset]).push_back(c[i]);
		}
		else {
			std::get<0>(tree[offset]).push_back(c[i]);
			dendrogram_to_tree(dendrogram, c[i], tree);
		}	
	}
}


hdl_tree_t convert_dendrogram_to_tree(const dendrogram_t dendrogram)
{
	hdl_tree_t tree;
	unsigned long int dendrogram_root_id = dendrogram.size()*2;
	dendrogram_to_tree(dendrogram, dendrogram_root_id, tree);
	cout << "DEBUG: " << dendrogram.size() << endl;
	return tree;
}




uint32_t add_num_of_leafs_to_dendrogram(dendrogram_t &dendrogram, unsigned long int merge)
{
	uint32_t num_leafs = dendrogram.size()+1;
	uint32_t merge_index = merge-num_leafs;

	if(merge < num_leafs)
		return 1;

	uint32_t n =
		add_num_of_leafs_to_dendrogram( dendrogram, std::get<0>(dendrogram[merge_index]) ) +
		add_num_of_leafs_to_dendrogram( dendrogram, std::get<1>(dendrogram[merge_index]) );
	std::get<2>(dendrogram[merge_index]) = n;
	return n;
}


dendrogram_t improve_igraph_clustering_merges(igraph_dendrogram_merges_t igraph_dendrogram)
{
	uint32_t dendrogram_root = igraph_dendrogram.size()*2;
	dendrogram_t dendrogram;

	for(unsigned long int i = 0; i<igraph_dendrogram.size(); i++) {
		dendrogram.push_back( std::make_tuple( get<0>(igraph_dendrogram[i]), get<1>(igraph_dendrogram[i]), -1 ) );
	}

	add_num_of_leafs_to_dendrogram(dendrogram, dendrogram_root);
	return dendrogram;
}



//vector<uint32_t> get_distance_vector(hdl_tree_t tree, uint32_t index



double cophenetic_correlation_v2(hdl_hier_t hdl, hdl_tree_t tree, int request_nthreads)
{

	cout << "DEBUG: " << hdl.size() << endl;
	return 111;

	int num_of_cells = 0;
	// Print tree
	for(hdl_tree_t::const_iterator it = tree.begin(); it != tree.end(); ++it)
	{
		const vector<uint32_t> *sons = &(it->second.first);
		const vector<uint32_t> *cells = &(it->second.second);

		cout << it->first << ((sons->size() == 0)? " - NODE\n" : " - LEAF\n");

		//assert (sons->size() == 0 || cells->size() == 0);
		for(long unsigned int i=0; i<sons->size(); i++) cout << (*sons)[i] << ", ";
		for(long unsigned int i=0; i<cells->size(); i++) cout << (*cells)[i] << ", ";
		cout << endl;

		num_of_cells += cells->size();
	}
	cout << "@@@@@@@@@@@@@@@@@@ " << num_of_cells << endl;




	double tree_distance_mean = 0;
	uint64_t tree_vec_size = 0;


	for(hdl_tree_t::const_iterator it = tree.begin(); it != tree.end(); ++it)
	{
		const vector<uint32_t> *sons = &(it->second.first);
		const vector<uint32_t> *cells = &(it->second.second);

		cout << it->first << ((sons->size() == 0)? " - NODE\n" : " - LEAF\n");

		//assert (sons->size() == 0 || cells->size() == 0);
		for(long unsigned int i=0; i<sons->size(); i++) cout << (*sons)[i] << ", ";
		for(long unsigned int i=0; i<cells->size(); i++) cout << (*cells)[i] << ", ";
		cout << endl;

		num_of_cells += cells->size();
	}

/*

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_tree_distance_mean = 0;
		uint64_t l_tree_vec_size = 0;


typedef map <int, pair< vector<uint32_t>, vector<uint32_t>>> hdl_tree_t;

		//for(long unsigned int i=id; i<_tree.size(); i+=nthrds) {
		for(long unsigned int i=id; i<tree.size(); i+=nthrds) {
			vector<uint32_t> vec = tree.get_distance_vector(i);
			for (long unsigned int j = 0; j < vec.size(); j++) {
				l_tree_distance_mean += vec[j];
				++l_tree_vec_size;
			}
			if( (i & 0xFFF)  == 0)
				cout << id << ": " << i << "/" << _tree.size() << endl;
		}
		#pragma omp critical
		{
			tree_distance_mean += l_tree_distance_mean;
			tree_vec_size += l_tree_vec_size;
		}

	}

	cout << "#threads: " <<  nthreads <<  endl;
	tree_distance_mean /= tree_vec_size;
	cout << "tree_distance_mean: " <<  tree_distance_mean <<  endl;


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

		for (long unsigned int i = id; i < hdl.size(); i+=nthrds) {
			for (long unsigned int j = i+1; j < hdl.size(); j++) {
				l_hdl_hierarchy_distance_mean += hdl_hierarchy_distance(hdl[i].c_str(), hdl[j].c_str());
				++l_hdl_vec_size;
			}
			if( (i & 0xFFF)  == 0)
				cout << id << ": " << i << "/" << hdl.size() << endl;
		}

		#pragma omp critical
		{
			hdl_hierarchy_distance_mean += l_hdl_hierarchy_distance_mean;
			hdl_vec_size += l_hdl_vec_size;

		}

	}

	cout << "#threads: " <<  nthreads <<  endl;
	hdl_hierarchy_distance_mean /= hdl_vec_size;
	cout << "hdl_hierarchy_distance_mean: " <<  hdl_hierarchy_distance_mean <<  endl;


	if(dendrogram_vec_size != hdl_vec_size) {
		cout << "ERROR: " << dendrogram_vec_size << " != " << hdl_vec_size << endl;
		return -1;
	}



	double c_numerator, c_denominator1, c_denominator2;

	omp_set_num_threads(request_nthreads);
	#pragma omp parallel
	{

		int id = omp_get_thread_num();
		int nthrds = omp_get_num_threads();
		if(id == 0) nthreads = nthrds;

		double l_c_numerator = 0, l_c_denominator1 = 0, l_c_denominator2 = 0;

		for(long unsigned int i=id; i<_dendrogram.size(); i+=nthreads) {
			vector<uint32_t> dendrogram_vec = dendrogram.get_distance_vector(i);
			vector<uint32_t> hdl_vec;
			for (long unsigned int j = i+1; j < hdl.size(); j++)
				hdl_vec.push_back(hdl_hierarchy_distance(hdl[i].c_str(), hdl[j].c_str()));

			for(long unsigned int j=0; j<dendrogram_vec.size(); j++) {
				double hdl_tmp = hdl_vec[j]-hdl_hierarchy_distance_mean;
				double den_tmp = dendrogram_vec[j]-dendrogram_distance_mean;

				l_c_numerator += hdl_tmp*den_tmp;
				l_c_denominator1 += pow(hdl_tmp,2);
				l_c_denominator2 += pow(den_tmp,2);
			}
			if( (i & 0xFFF)  == 0)
				cout << id << ": " << i << "/" << _dendrogram.size() << endl;
		}

		#pragma omp critical
		{
			c_numerator += l_c_numerator;
			c_denominator1 += l_c_denominator1;
			c_denominator2 += l_c_denominator2;
		}
	}
	cout << "#threads: " <<  nthreads <<  endl;

	double cophenetic_factor = c_numerator / sqrt(c_denominator1*c_denominator2);
	cout << "RESULT: " << cophenetic_factor << endl;
	return cophenetic_factor;
	*/
	return 0;
}











PYBIND11_MODULE(cophenetic_module, handle) {
	handle.doc() = "cophenetic_correlation...";
	handle.def("cophenetic_correlation", &cophenetic_correlation);

	handle.def("cophenetic_correlation_v2", &cophenetic_correlation_v2);
	handle.def("improve_igraph_clustering_merges", &improve_igraph_clustering_merges);
	handle.def("convert_dendrogram_to_tree", &convert_dendrogram_to_tree);
}


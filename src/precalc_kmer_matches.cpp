#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>
#include <string>
#include <list>
#include <utility>
#include <boost/functional/hash.hpp>
#include <fstream>


using namespace sdsl;

//what to do with Ns?

//void generate_all_kmers(std::vector<uint8_t> letters, std::vector<uint8_t>& substr, int k, int n, std::vector<std::vector<uint8_t>>& kmers);
//void get_kmers(char *kmerfile, std::vector<std::vector<uint8_t>>& kmers);

void precalc_kmer_matches (csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa, int k,   
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_idx, 
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_idx_rev,
			   SiteOverlapTracker* global_tracker_array,//prealloced array
			   SiteMarkerArray * prg_sites,
			   sequence_map<std::vector<uint8_t>, uint32_t>& kmer_to_tracker_index,
			   std::vector<int> mask_a, uint64_t maxx, 
			   sequence_set<std::vector<uint8_t>>& kmers_in_ref, char * kmerfile) 
{

	std::list<std::vector<std::pair<uint32_t, std::vector<int>>>> temp2;		  
	interval_list temp;
	bool first_del;

	ifstream kfile;
	string line;
	kfile.open(kmerfile);
	int kmer_counter=0;
	while (std::getline(kfile,line))
	{
		std::vector<uint8_t> kmer;
		for (auto c: line)
			switch (c)
			{
				case 'A': case 'a': kmer.push_back(1);break;
				case 'C': case 'c': kmer.push_back(2);break;
				case 'G': case 'g': kmer.push_back(3);break;
				case 'T': case 't': kmer.push_back(4);break;
			}

		kmer_idx[kmer]=temp;
		kmer_idx_rev[kmer]=temp;
		//kmer_sites[kmer]=temp2;
		kmer_to_tracker_index[kmer]=kmer_counter;

		first_del=false;
		std::vector<uint8_t>::iterator res_it=
		  bidir_search_bwd(csa,0,csa.size(),0,csa.size(),
				   kmer.begin(),kmer.end(),
				   kmer_idx[kmer],kmer_idx_rev[kmer],
				   &global_tracker_array[kmer_counter],
				   prg_sites,
				   mask_a,maxx,first_del);

		if (!first_del) kmers_in_ref.insert(kmer);
		kmer_counter++;
	}
}


//void generate_all_kmers(std::vector<uint8_t> letters, std::vector<uint8_t>& substr, int k, int n, std::vector<std::vector<uint8_t>>& kmers)
//{
//  if (k==1) {
//    for (int j = 0; j < n; j++)
//      substr.push_back(letters[j]);
//      kmers.push_back(substr);
//  }
//  else {
//    for (int jj = 0; jj < n; jj++)
//      substr.push_back(letters[jj]);
//      generate_all_kmers(letters, substr, k-1, n,kmers);
//  }
//}

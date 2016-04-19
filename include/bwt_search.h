#ifndef __BWT_SEARCH_H_INCLUDED__
#define __BWT_SEARCH_H_INCLUDED__

#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <boost/functional/hash.hpp> 
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <sitemarker.hpp>

using namespace sdsl;
using namespace std;

template < typename SEQUENCE > struct seq_hash
			       {
				 std::size_t operator() ( const SEQUENCE& seq ) const
				 {
				   std::size_t hash = 0 ;
				   boost::hash_range( hash, seq.begin(), seq.end() ) ;
				   return hash ;
				 }
};

template < typename SEQUENCE, typename T >

using sequence_map = std::unordered_map< SEQUENCE, T, seq_hash<SEQUENCE> > ;

template < typename SEQUENCE >

using sequence_set = std::unordered_set< SEQUENCE, seq_hash<SEQUENCE> > ;


typedef std::vector<std::pair<uint64_t,uint64_t>> interval_list;
typedef std::vector<SiteOverlapTracker> site_tracker_list;



csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa_constr(std::string fname, std::vector<std::vector<int>>& covgs, char* int_al_fname, char* memory_log_fname, char* csa_file, bool fwd);

void precalc_kmer_matches (csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa, int k,   
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_sa_interval, //proper (nonzero)
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_sa_interval_rev,
			   sequence_map<std::vector<uint8_t>, site_tracker_list>& kmer_tracker,
			   std::vector<int> mask_a, uint64_t maxx, 
			   sequence_set<std::vector<uint8_t>>& kmers_in_ref, char * kmerfile,
			   SiteInfo* si);
 
uint64_t parse_masks(std::vector<uint64_t>& mask_s, std::vector<int>& mask_a, std::string sites_fname, std::string alleles_fname, std::vector<std::vector<int>>& covgs);

uint64_t bidir_search(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa, 
		      uint64_t& left, uint64_t& right, 
		      uint64_t& left_rev, uint64_t& right_rev, 
		      uint8_t c);


void get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa,
		  uint64_t num_idx,
		  uint32_t num, bool last,
		  std::vector<int>& mask_a,
		  std::vector<SiteOverlapTracker>::iterator& it_s,
		  SiteInfo* site_info);

bool skip(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa,
                      uint64_t& left, uint64_t& right,
                      uint64_t& left_rev, uint64_t& right_rev,
	              uint32_t num);

std::vector<uint8_t>::iterator bidir_search_bwd(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa,
						uint64_t left, uint64_t right,
						uint64_t left_rev, uint64_t right_rev,
						std::vector<uint8_t>::iterator pat_begin, 
						std::vector<uint8_t>::iterator pat_end,
						interval_list& sa_intervals, 
						interval_list& sa_intervals_rev,
						interval_list& sa_intervals_temp, 
						interval_list& sa_intervals_rev_temp,
						site_tracker_list& site_trackers,//one per proper interval, pre-reserved
						site_tracker_list& site_trackers_temp,//one for all intervals during search
						std::vector<int> mask_a, uint64_t maxx, bool& first_del,
						SiteInfo* site_info
						);

std::vector<uint8_t>::iterator bidir_search_fwd(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa,
						uint64_t left, uint64_t right,
						uint64_t left_rev, uint64_t right_rev,
						std::vector<uint8_t>::iterator pat_begin, std::vector<uint8_t>::iterator pat_end,
						std::list<std::pair<uint64_t,uint64_t>>& sa_intervals,
						std::list<std::pair<uint64_t,uint64_t>>& sa_intervals_rev,
						std::list<std::vector<std::pair<uint32_t, std::vector<int>>>>& sites,
						std::vector<int> mask_a, uint64_t maxx, bool& first_del);

#endif

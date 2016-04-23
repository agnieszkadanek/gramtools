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
#include "sitemarker.hpp"

using namespace sdsl;

//what to do with Ns?

void precalc_kmer_matches (csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2>& csa, int k,   
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_sa_interval, //proper (nonzero)
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_sa_interval_rev,
			   sequence_map<std::vector<uint8_t>, site_tracker_list>& kmer_tracker,
			   std::vector<int>& mask_a, uint64_t maxx, 
			   sequence_set<std::vector<uint8_t>>& kmers_in_ref, char * kmerfile,
			   SiteInfo* si) 
{
  
  //set up some temporary lists for use during this function
  interval_list temp, temp_rev;
  temp.reserve(10000);
  temp_rev.reserve(10000);
  site_tracker_list site_trackers_temp;
  site_trackers_temp.reserve(10000);

  site_tracker_list::iterator it_s;
  it_s = site_trackers_temp.begin();
  

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
      kmer_sa_interval[kmer].reserve(10000);
      kmer_sa_interval_rev[kmer].reserve(10000);
      
      first_del=false;
      std::vector<uint8_t>::iterator res_it=
	bidir_search_bwd(csa,0,csa.size(),0,csa.size(),
			 kmer.begin(),kmer.end(),
			 kmer_sa_interval[kmer], kmer_sa_interval_rev[kmer],
			 temp, temp_rev,
			 kmer_tracker[kmer],
			 site_trackers_temp,
			 mask_a,maxx,first_del, si);


      if (!first_del) kmers_in_ref.insert(kmer);
      ++it_s;
      temp.resize(0);
      temp_rev.resize(0);
    }
}


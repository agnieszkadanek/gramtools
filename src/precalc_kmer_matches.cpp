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

void precalc_kmer_matches (csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa, int k,   
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_sa_interval, //proper (nonzero)
			   sequence_map<std::vector<uint8_t>, interval_list>& kmer_sa_interval_rev,
			   std::vector<SiteOverlapTracker>& site_trackers,//one per proper (nonempty) interval
			   std::vector<SiteOverlapTracker>& site_trackers_temp,//one per interval
			   sequence_map<std::vector<uint8_t>, uint32_t>& kmer_index,//which-th kmer is it
			   std::vector<int> mask_a, uint64_t maxx, 
			   sequence_set<std::vector<uint8_t>>& kmers_in_ref, char * kmerfile) 
{
  
  std::vector<SiteOverlapTracker>::iterator it_s;
  //site_trackers_temp.clear();
  it_s = site_trackers_temp.begin();
  
  interval_list temp, temp_rev;
  temp.reserve(10000);
  temp_rev.reserve(10000);

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
      
      kmer_index[kmer]=kmer_counter;
      
      first_del=false;
      std::vector<uint8_t>::iterator res_it=
	bidir_search_bwd(csa,0,csa.size(),0,csa.size(),
			 kmer.begin(),kmer.end(),
			 kmer_sa_interval[kmer], kmer_sa_interval_rev[kmer]
			 temp, temp_rev,
			 site_trackers_temp,
			 mask_a,maxx,first_del);


      if (!first_del) kmers_in_ref.insert(kmer);
      kmer_counter++;
      ++it_s;
      temp.clear();
      temp_rev.clear();
    }
}


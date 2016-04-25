#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
//#include "bwt_search.h"
#include <tuple>
#include <cstdint>

using namespace sdsl;
using namespace std;

uint64_t bidir_search_with_ranks(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2>& csa, 
		      uint64_t& left, uint64_t& right, 
		      uint64_t& left_rev, uint64_t& right_rev, 
		      uint8_t c, unordered_map<uint8_t,vector<uint64_t>>& rank_all)
{
  assert(left < right); 
  assert(right <= csa.size());
  assert((c>0) & (c<5));  

  uint64_t c_begin = csa.C[csa.char2comp[c]];

  left=c_begin+rank_all[c][left-1];
  right=c_begin+rank_all[c][right-1];
  assert(right>=left);

  //need to calc s and b from rank_all
  //left_rev  = left_rev + s;
  //right_rev = right_rev - b + 1;
  //assert(right_rev-left_rev == right-left);

  return right-left;
}

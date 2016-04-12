#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>
#include "sitemarker.hpp"

using namespace sdsl;

void get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa,
		  uint64_t num_idx,
		  uint32_t num, bool last,
		  std::vector<int>& allele,
		  std::vector<int> mask_a,
		  SiteOverlapTracker* tracker)
{
  uint32_t site;
  
  if (num%2==1) {
    site=(num-5)/2;;
    if (!last)
      {
	//allele.push_back(1);
	tracker.push(site, 1 );
      }
  }
  else {
    site=((num-1)-5)/2;
    tracker.push(site, mask_a[csa[num_idx]] );
    //allele.push_back(mask_a[csa[num_idx]]);
  }
  
}  

//removed left and right args because they only form a proper interval is when last is true; but then they can't help with location because allele is unknown and we're just addin site

  //problems:
  //avoid mask_s
  //when have matches to the right of both odd numbers

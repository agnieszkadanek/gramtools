#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include "bwt_search.h"
#include <tuple>
#include <cstdint>
#include "sitemarker.hpp"

using namespace sdsl;

void get_location(csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2>& csa,
		  uint64_t num_idx,
		  uint32_t num, bool last,
		  std::vector<int>& mask_a,
		  std::vector<SiteOverlapTracker>::iterator& it_s,
		  SiteInfo* site_info)

{

  uint32_t site;
  uint32_t allele;
  if (num%2==1) 
    {
      site=(num-5)/2;;
      if (!last)
	{
	  allele=0;
	}
      else
	{
	  return;
	}
    }
  else 
    {
      site=((num-1)-5)/2;
      allele=mask_a[csa[num_idx]];
      printf("got allele %d\n", allele);
    }
  
  (*it_s).push(site, allele, site_info);

}  

//removed left and right args because they only form a proper interval is when last is true; but then they can't help with location because allele is unknown and we're just addin site

  //problems:
  //avoid mask_s
  //when have matches to the right of both odd numbers

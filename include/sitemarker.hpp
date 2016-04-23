#ifndef __SITEMARKER__
#define __SITEMARKER__

#include <vector>
#include <boost/dynamic_bitset.hpp>


//this object gets allocated once - parse a file containing one line per alpha
// the sites vector has alphabet+1 entries. at index 0,1,2,3,4 - just have NULL pointers or something
// then you index each site with the number it has in the linPRG. 
class SiteInfo
{
public:
  SiteInfo(std::string sitefile);
  ~SiteInfo();
  
  //suppose you want to get a pointer to the SiteMarker for site 56
  //and at the same time set alleles 1,5,6 to 1 (meaning this read overlaps alleles 1,5,6:
  //note this will zero all other bits.
  uint32_t get_num_alleles(uint32_t site_id);
  uint32_t get_num_sites();
private:
  std::vector<uint32_t> allele_counts;
  uint32_t num_sites;//literally the number of sites, nothing to do with int-alphabet/even/odd stuff
};


//what is the thing we use when tracking which sites/alleles a read crosses?
//we will have one of these per interval
class SiteOverlapTracker
{
public:
  SiteOverlapTracker(): alleles(2, boost::dynamic_bitset<>(2))
  {
    //sites.reserve(2);
  }
  SiteOverlapTracker(const SiteOverlapTracker&);
  void push(uint32_t site_id, uint32_t allele, SiteInfo* si);
  void set(uint32_t site_id, uint32_t allele, SiteInfo* si);
  void clear();
  std::vector<uint32_t> sites;
  std::vector<boost::dynamic_bitset<>> alleles;  
};



#endif


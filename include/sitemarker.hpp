#ifndef __SITEMARKER__
#define __SITEMARKER__

#include <vector>
#include <boost/dynamic_bitset.hpp>


/*
class SiteMarker
{
public:
  int int_alphabet_odd_id; //odd number corresponding to this site
  int site_index;//which-th site is this. Runs from 0 to num_sites-1
  

  //site index is like 0,1,2,3,....
  //odd_id is like 5,7,9,11....
  //odd_id = 2*site_index+5
  SiteMarker(int odd_id, int num_alleles): int_alphabet_odd_id(odd_id), 
					   site_index((odd_id-5)/2), 
					   alleles(num_alleles) {};

  void zero_all_alleles();
  void set_allele(uint32_t i);//set it to 1 if we cross this allele
  void set_these_alleles(std::vector<int> v); //to set 5th allele, do    alleles[5]=1;
  int count_set_alleles();
  int get_num_alleles(); // this returns alleles.size();
  int get_allele_bit(uint32_t i);
  void print_all_info();

private:



};

*/

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
  int num_sites;//literally the number of sites, nothing to do with int-alphabet/even/odd stuff
};



  




//what is the thing we use when tracking which sites/alleles a read crosses?
//we will have one of these per interval
class SiteOverlapTracker
{
public:
  SiteOverlapTracker(): valid(true),
			alleles(2, boost::dynamic_bitset<>(2));//by default, 2 sites each with 2 alleles
  SiteOverlapTracker(const SiteOverlapTracker&);
  void push(int site_id, int allele, SiteInfo* si);
  void clear();
  boolean is_valid();
  void invalidate();
private:
  boolean valid;
  std::vector<uint32_t> sites;
  std::vector<boost::dynamic_bitset<>> alleles;  
};


class SiteOverlapTrackerArray
{
public:
  std::vector<SiteOverlapTracker*> site_trackers;
  uint32_t get_next(uint32_t current_inde);
};


#endif


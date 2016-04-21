/*
 * 
 */
#include <vector>
#include <string>
#include <cstdlib>
#include "sitemarker.hpp"
#include <algorithm>
#include <fstream>


//pass in a file which has one column
//the i-th row = number of alleles in i-th site.
SiteInfo::SiteInfo(std::string sitefile)
{
  //count lines in file so know how much to alloc
  std::ifstream inFile(sitefile); 
  num_sites = std::count(std::istreambuf_iterator<char>(inFile), 
			 std::istreambuf_iterator<char>(), '\n');
  inFile.close();
  printf("Found %d sites\n", num_sites);

  //collect info on how many alleles in each site from input file
  allele_counts.reserve(num_sites);

  std::ifstream fs(sitefile);
  int i=0;
  for(std::string line; std::getline(fs, line); )
    {
      allele_counts.push_back(stoi(line));
      i++;
    }
}

SiteInfo::~SiteInfo()
{
}

uint32_t SiteInfo::get_num_alleles(uint32_t site_id)
{
  if (site_id>num_sites-1)
    {
      printf("Calling site id %d which is too big. Expected to be < %d. Also - need to put proper exception handling in\n",
	     site_id, num_sites-1);
      exit(1);
    }
  return allele_counts[site_id];
}

uint32_t SiteInfo::get_num_sites()
{
  return num_sites;
}


SiteOverlapTracker::SiteOverlapTracker(const SiteOverlapTracker& from):sites(from.sites), alleles(from.alleles)					 
{
}



void SiteOverlapTracker::push(uint32_t site_id, uint32_t allele, SiteInfo* si)
{
  if ((sites.size()>0) && (site_id==sites.back()))
    {
      alleles.back()[allele]=1;//flick the bit for this allele
    }
  else//either no sites at all yet or a new site
    {
      printf("Push %d\n", site_id);
      sites.push_back(site_id);
      uint32_t num = si->get_num_alleles(site_id);
      if (allele>num)
	{
	  printf("Trying to modify an allele %d with out of bounds index beyond %d\n", allele, num-1);
	  exit(1);
	}
      else if (allele<0)
	{
	  printf("Trying to modify negative allele %d\n", allele);
	  exit(1);
	}
      alleles.push_back(boost::dynamic_bitset<>(num) );
      // TEMPORARY ZAM alleles.back()[allele]=1;
    }
}


void SiteOverlapTracker::clear()
{
  sites.clear();
  alleles.clear();
}


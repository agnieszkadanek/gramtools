#include "sdsl/suffix_arrays.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <fstream>
#include <iostream>
#include <cstdint>
#include <time.h> 
#include <vector>  
#include "bwt_search.h"
#include <seqread.hpp>

using namespace std;  
using namespace sdsl;

void timestamp(); 

//argv[1] - file containing linear prg
//argv[2] - file where CSA is stored 
//argv[3] - file containing reads to be mapped (one read per line)
//argv[4] - file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which site
//argv[5] - file containing mask over the linear prg that indicates at each position whether you are inside a site and if so, which allele
//argv[6] - name of output file where coverages on each allele are printed
//argv[7] - name of output file where reads that have been processed are printed
//argv[8] - name of binary file where the prg in integer alphabet is written
//argv[9] - memory log file for CSA
//argv[10] - size of precalculated kmers
//argv[11] - kmer file
//argv[12] - file, one line per site, each line just has an integer - number alleles at that site
//argv[13] - number of kmers in kmer file
int main(int argc, char* argv[]) {
  
  std::vector<uint64_t> mask_s;
  std::vector<int> mask_a;
  std::vector<std::vector<int> > covgs;
  std::string q;
  
  SeqRead inputReads(argv[3]); 
  std::ofstream out(argv[6]); 
  std::ofstream out2(argv[7]);
  
  SiteInfo* si = new SiteInfo(argv[12]);

  //associative array, key->value are kmer->list of BWT intervals
  sequence_map<std::vector<uint8_t>, interval_list> kmer_sa_intervals,kmer_sa_intervals_rev;
  kmer_sa_intervals.reserve(atoi(argv[13]));
  kmer_sa_intervals_rev.reserve(atoi(argv[13]));
  
  //assoc array, key-> value are kmer->list of site_overlap_trackers
  sequence_map<std::vector<uint8_t>, site_tracker_list> kmer_to_trackerlist;
  kmer_to_trackerlist.reserve(atoi(argv[13]));
  
  sequence_set<std::vector<uint8_t>> kmers_in_ref;
  
  //not using mask_s anymore?
  uint64_t maxx=parse_masks(mask_s,mask_a,argv[4],argv[5],covgs);
  
  timestamp();
  cout<<"CSA construction"<<endl;
  csa_wt<wt_int<bit_vector,rank_support_v5<>>,2,2> csa=csa_constr(argv[1],covgs,argv[8],argv[9],argv[2],true);
  timestamp();
  
  std::vector<std::vector<string> > site_reads(covgs.size(),std::vector<string>(1));
  int no_reads=0;
  uint64_t no_mapped=0;
  int inc=0;
  int no_occ=0;
  bool first_del=false;
  
  int k=atoi(argv[10]); //verify input
  precalc_kmer_matches(csa,k,
		       kmer_sa_intervals,kmer_sa_intervals_rev,			     
		       kmer_to_trackerlist,
		       mask_a,maxx,kmers_in_ref,argv[11], si);
  timestamp();
  
  //setup variables for main bidirectional search
  interval_list sa_intervals, sa_intervals_rev, temp, temp_rev;
  site_tracker_list site_trackers, site_trackers_temp;
  sa_intervals.reserve(10000);
  sa_intervals_rev.reserve(10000);
  temp.reserve(10000);
  temp_rev.reserve(10000);
  site_trackers.reserve(10000);


  interval_list::iterator it, it_rev, it_temp;
  std::vector<uint8_t>::iterator res_it;
  
  for (auto q: inputReads)
    {
      // If you declare p inside the while scope, it is destroyed/created automatically in every loop
      std::vector<uint8_t> p;
      
      //logging
      if (!(inc++%10)) { out2<<no_reads<<endl; }
      
      //add N's
      
      for (int i=0,seqlen=strlen(q->seq);i<seqlen;i++) {
	if (q->seq[i]=='A' or q->seq[i]=='a') p.push_back(1);
	if (q->seq[i]=='C' or q->seq[i]=='c') p.push_back(2);
	if (q->seq[i]=='G' or q->seq[i]=='g') p.push_back(3);
	if (q->seq[i]=='T' or q->seq[i]=='t') p.push_back(4);
      }
      
      std::vector<uint8_t> kmer(p.begin()+p.size()-k,p.end()); //is there a way to avoid making this copy?
      
      //do we really need this .find()? It's going to be (worst case) linear in number of kmers.
      if (kmer_sa_intervals.find(kmer)!=kmer_sa_intervals.end() && 
	  kmer_sa_intervals_rev.find(kmer)!=kmer_sa_intervals_rev.end() && 
	  kmer_to_trackerlist.find(kmer)!=kmer_to_trackerlist.end()) 
	{
	  sa_intervals=kmer_sa_intervals[kmer];
	  sa_intervals_rev=kmer_sa_intervals_rev[kmer];
	  site_trackers=kmer_to_trackerlist[kmer];	
	  
	  it=sa_intervals.begin();
	  it_rev=sa_intervals_rev.begin();
	  
	  if (kmers_in_ref.find(kmer)!=kmers_in_ref.end()) first_del=false;
	  else first_del=true;
	  temp.clear();
	  temp_rev.clear();
	  res_it=bidir_search_bwd(csa, (*it).first, (*it).second, 
				  (*it_rev).first, (*it_rev).second, 
				  p.begin(),p.begin()+p.size()-k, 
				  sa_intervals, sa_intervals_rev, 
				  temp, temp_rev,
				  site_trackers, site_trackers_temp,
				  mask_a, maxx, first_del, si);
	  
	  no_occ=0;
	  for (it=sa_intervals.begin();it!=sa_intervals.end();++it)
	    no_occ+=(*it).second-(*it).first;
	  
	  sa_intervals.clear();
	  sa_intervals_rev.clear();
	  temp.clear();
	  temp_rev.clear();
	  site_trackers.clear();
	  site_trackers_temp.clear();
	}
      else 
	{
	  no_occ=0;
	}
      
      out<<no_occ<<" ";
      //clear p, sa_intervals etc
      p.clear();
      
      no_reads++;
    }
  timestamp();
  
  out.close();
  out2.close();
  return(0);
}

void timestamp(){
	time_t ltime;
	ltime = time(NULL);
	printf("\n-----\n%s",asctime(localtime(&ltime)));
	fflush(stdout);
}

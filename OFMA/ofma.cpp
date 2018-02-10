#include <glob.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <map>
#include <math.h>
using namespace std;


//#define CDEBUG 1

// util for debug
#ifdef CDEBUG
#define LOG(m,x,y,z) log_message(m,x,y,z)
#else
#define LOG(m,x,y,z)
#endif

void log_message(string m, double x, double y, double z) {
  cerr << m << "," << x << "," << y  << "," << z << "," << "\n";
}


// main
int main (int argc, char* argv[])
{
  // parameters
  if(argc != 3) {
    return 1;
  }
  const char* path = argv[1];
  int64_t max_size(atoll(argv[2]));

  ifstream infile;
  infile.open(path);

  uint64_t t, id, size;

  std::string line;

  unordered_map<uint64_t,uint64_t> cache;
  int64_t cache_size = 0;
  
  uint64_t reqs = 0, hits = 0;
  uint64_t uint64_t rbytes = 0, hbytes = 0;

  bool logStatistics=true; // start empty, no warmup

  uint64_t counter = 0;

  while (infile >> t >> id >> size)
    {
      if (logStatistics)
	{
	  reqs++;
	  rbytes += size;
	}
      if (cache.count(id) > 0)
	{
	  // hit: already in cache
	  if (logStatistics)
	    {
	      hits ++;
	      hbytes += size;
	    }
	  LOG("hit",id,size,cache_size);
	}
      else
	{
	  // miss: admit and follow Belady's rule for eviction
	  cache[id] = size;
	  cache_size += size;
	  LOG("miss",id,size,cache_size);
	  if (cache_size > max_size) {
	    // find next requests for all
	    map<uint64_t,uint64_t> futureReqs;
	    unordered_map<uint64_t,uint64_t> tobefound (cache); // copy cache map
	    // save current position
	    int len = infile.tellg();
	    uint64_t reqdist = 0;
	    while (!tobefound.empty())
	      {
		if (++reqdist > 100000)
		  {
		    for (auto it: tobefound)
		      {
			futureReqs[reqdist] = get<0>(it);
			// bug?:			futureReqs[reqdist++] = get<0>(it);
			//			LOG("forloop",reqdist-1,get<0>(it),get<1>(it));
		      }
		    break;
		  }

		// peak ahead and if found cached object: record how many reqdist
		infile >> t >> id >> size;
		unordered_map<uint64_t,uint64_t>::const_iterator got = tobefound.find (id);
		if ( got != tobefound.end() )
		  {
		    //		    LOG("found",id,tobefound[id],cache[id]);
		    futureReqs[reqdist] = id;
		    tobefound.erase(got);
		  }
	      }
	    // remove objects with requests furthest in future (highest reqdist) from cache
	    map<uint64_t,uint64_t>::const_iterator it = futureReqs.end();
	    unordered_map<int,int> classlist;
	    while (--it != futureReqs.begin())
	      {
		uint64_t eid = (it)->second;
		uint64_t esize = cache[eid];
		// evict up to two object from each size class
		if ( ++classlist[floor(log2(esize))] < 3)
		  {
		    LOG("evict",floor(log2(esize)),eid,esize);
		    cache_size -= esize;
		    cache.erase(eid);
		  }
		futureReqs.erase(it);
	      }
	    LOG("cachesize",max_size,cache_size,double(cache_size)/max_size);
	    // return to saved position
	    infile.seekg(len ,std::ios_base::beg);
	  }
	}
      if(counter++ > 1000) {
	cout << "OFMA " << max_size << " hitc " << hits << " reqc " << reqs << " OHR " << double(hits)/reqs << "\n";
	counter=0;
      }
    }
  infile.close();

  cout << "OFMA " << max_size << " hitc " << hits << " reqc " << reqs << " OHR " << double(hits)/reqs << "\n";
  
  return 0;
}

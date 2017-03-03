#include <glob.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <math.h>
using namespace std;

typedef tuple<long, long, long> tuple_t; // objectid, lastr?, size
typedef vector<tuple_t>::iterator vector_iterator_t;

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

// original belady with search boundary



// main
int main (int argc, char* argv[])
{
  // parameters
  if(argc != 3) {
    return 1;
  }
  const char* path = argv[1];
  uint64_t max_size(atoll(argv[2]));

  const long lineCount = 1000000;
  long currentLine = 0;

  ifstream infile;
  infile.open(path);

  long t, id, size;

  // store whether object is requested again
  vector<tuple_t> req_vec;
  req_vec.reserve(lineCount+1);
  unordered_map<long, vector_iterator_t> prev_acc;
  unordered_map<long, long> size_map;

  cerr << "running...\n";

  while (infile >> t >> id >> size)
    {
      req_vec.push_back (tuple_t(id,0,size));
      // check if need to correct previous entry
      if(prev_acc.count(id)>0)
	{
	  get<1>(*(prev_acc[id])) = currentLine; // store when future request hapens
	}
      // set reference to last added element
      prev_acc[id] = (--req_vec.end());
      //
      // progress bar
      if(++currentLine % (lineCount/25) == 0)
	cerr << double(currentLine)/lineCount << '\n';
      if (currentLine > lineCount)
	break;
    }
  infile.close();

  cerr << "----------------\n";

  long reqs = 0, hits = 0;
  long long rbytes = 0, hbytes = 0;

  unordered_map<long,long> cache;
  map<long,long> cache_future;
  unordered_map<long,long> reverse_cache_future;
  long long cache_size = 0;

  currentLine = 0;
  bool logStatistics=false;

  for (vector_iterator_t it = req_vec.begin() ; it != req_vec.end(); it++)
    {
      // start statistics after warm up
      if (!logStatistics && t > 0) //two days warm up
	{
	  cerr << "statistics started.\n";
	  logStatistics = true;
	}
      if (logStatistics)
	{
	  reqs++;
	  rbytes += size;
	}
      id = get<0>(*it);
      long future_req = get<1>(*it);
      size = get<2>(*it);
      if (cache.count(id) > 0 && size == cache[id])
	{
	  // hit and consistent object size
	  if (logStatistics)
	    {
	      hits ++;
	      hbytes += size;
	    }
	}
      else
	{
	  // miss
	  // check if need to correct inconsistent object size
	  if (size != cache[id])
	    {
	      cache_size -= cache[id];
	      cache.erase(id);
	      cache_future.erase(reverse_cache_future[id]);
	      reverse_cache_future.erase(id);
	    }
	  // only admit if there's a next request for same object in the trace
	  if (future_req > 0)
	    {
	      cache_future[future_req] = id; // map future_req -> id
	      reverse_cache_future[id] = future_req; // map id -> future_req
	      cache[id] = size; // map id -> size
	      cache_size += size;
	      while (cache_size > max_size)
		{
		  // evict objects with requests furthest in future (highest cache_future)
		  map<long,long>::const_iterator it = cache_future.end();
		  long eid = (--it)->second; // id of object
		  long esize = cache[eid]; // size of object
		  cache_size -= esize;
		  cache.erase(eid);
		  cache_future.erase(it);
		  reverse_cache_future.erase(eid);
		}
	    }
	}
      //
      // progress bar
      if(++currentLine % (lineCount/25) == 0)
	cerr << double(currentLine)/lineCount << '\n';
    }

  cout << "Belady " << max_size << " hitc " << hits << " reqc " << reqs << " OHR " << double(hits)/reqs << "\n";
  
  return 0;
}

#include <iostream>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <list>
#include <tuple>
#include <assert.h>
#include <random>


using namespace std;

typedef tuple<long, long> object_t; // objectid, size
typedef list<object_t>::iterator list_iterator_t;
typedef multimap<long double, object_t>::iterator map_iterator_t;

// set to enable cache debugging
//#define CDEBUG 1
//#define TDEBUG 1

// util for debug
#ifdef CDEBUG
#define LOG(m,x,y,z) log_message(m,x,y,z)
#else
#define LOG(m,x,y,z)
#endif

void log_message(string m, double x, double y, double z) {
  cerr << m << "," << x << "," << y  << "," << z << "\n";
}

/*
  general cache class
*/

class Cache {
protected:
  // caching definitions
  const long long cache_size;
  long long current_size; // size of the cache
  // hit statistics
  long hits;
  long long bytehits;

  bool logStatistics;
  virtual void hit(long size) {
    if(logStatistics) {
      hits++;
      bytehits+=size;
    }
  }

public:
  Cache(long long cs) : cache_size(cs), current_size(0), hits(0), bytehits(0), logStatistics(false) {
  }

  virtual ~Cache(){};

  virtual bool request (const long cur_req, const long size) {return(false);}

  // statistics
  virtual void startStatistics() {
    logStatistics = true;
  }

  virtual void stopStatistics() {
    logStatistics = false;
  }

  virtual void resetStatistics() {
    hits = 0;
    bytehits = 0;
  }
  
  virtual long getHits() const {
    return(hits);
  }

  virtual bool getLogStatistics() const {
    return(logStatistics);
  }

  virtual long long getBytehits() const {
    return(bytehits);
  }

  virtual long long getCurrentSize() const {
    return(current_size);
  }

  virtual long long getCacheSize() const {
    return(cache_size);
  }
};

 
/*
  Static cache policy, which tracks
*/

class StaticCache: public Cache {
protected:
  unordered_set<long> cache_set;
  unordered_map<long, object_t> history_map;

public:
  StaticCache(long long cs): Cache(cs) {}

  ~StaticCache(){}

  bool request (const long cur_req, const long size) {
    // gather statistics
    if (history_map.count(cur_req) == 0)
      history_map[cur_req] = object_t(1,size);
    else
      (get<0>(history_map[cur_req]))++;
    // check whether in cache
    if(cache_set.count(cur_req)>0)
      {
	Cache::hit(size);
	return(true);
      }
    return(false);
  }

  void activateHistory () {
    //    cout << "history activated start\n";
    map<double, long> history_order;
    // sort request statistic by weight
    for (auto ait: history_map)
      history_order.insert(pair<double,long>(
					     double(get<0>(ait.second))/get<1>(ait.second), // weight: number of requests/ size
					     ait.first)); // id
    // clear old cache
    cache_set.clear();
    current_size = 0;
    cerr << "history size: " << history_order.size() << "\n";
    while (history_order.size()>0 && current_size <= cache_size) {
      map<double, long>::iterator highest_weight = (--history_order.end());
      const long highest_id = highest_weight->second;
      const long highest_size = get<1>(history_map[highest_id]);
      if(current_size + highest_size <= cache_size) {
	current_size += highest_size;
	cache_set.insert(highest_id);
      }
      history_order.erase(highest_weight);
    }
    history_map.clear();
    cerr << "history activated - cache size: " << double(current_size)/cache_size << "\n";
  }
};




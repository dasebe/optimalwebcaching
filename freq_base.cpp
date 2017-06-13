#include <iostream>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <list>
#include <tuple>
#include <cassert>
#include <random>


using namespace std;

typedef tuple<uint64_t, uint64_t> object_t; // objectid, size
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
    const uint64_t cache_size;
    uint64_t current_size; // size of the cache
    // hit statistics
    uint64_t hits;
    uint64_t bytehits;

    bool logStatistics;
    virtual void hit(uint64_t size) {
        if(logStatistics) {
            hits++;
            bytehits+=size;
        }
    }

public:
    Cache(uint64_t cs) : cache_size(cs), current_size(0), hits(0), bytehits(0), logStatistics(false) {
    }

    virtual ~Cache(){};

    virtual bool request (const uint64_t cur_req, const uint64_t size) {return(false);}

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
  
    virtual uint64_t getHits() const {
        return(hits);
    }

    virtual bool getLogStatistics() const {
        return(logStatistics);
    }

    virtual uint64_t getBytehits() const {
        return(bytehits);
    }

    virtual uint64_t getCurrentSize() const {
        return(current_size);
    }

    virtual uint64_t getCacheSize() const {
        return(cache_size);
    }
};

 
/*
  Static cache policy, which tracks
*/

class StaticCache: public Cache {
protected:
    unordered_set<uint64_t> cache_set;
    unordered_map<uint64_t, object_t> history_map;

public:
    StaticCache(uint64_t cs): Cache(cs) {}

    ~StaticCache(){}

    bool request (const uint64_t cur_req, const uint64_t size) {
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

    bool antiRequest (const uint64_t cur_req, const uint64_t size) {
        //        cerr << "anti " << cur_req << " " << size << "\n";
        //        assert(history_map.count(cur_req) > 0);
        assert(get<0>(history_map[cur_req]) > 0);
        (get<0>(history_map[cur_req]))--;
    }

    void activateHistory () {
        //    cout << "history activated start\n";
        map<double, uint64_t> history_order;
        // sort request statistic by weight
        for (auto ait: history_map)
            history_order.insert(pair<double,uint64_t>(
                                                   double(get<0>(ait.second))/get<1>(ait.second), // weight: number of requests/ size
                                                   ait.first)); // id
        // clear old cache
        cache_set.clear();
        current_size = 0;
        cerr << "history size: " << history_order.size() << "\n";
        while (history_order.size()>0 && current_size <= cache_size) {
            map<double, uint64_t>::iterator highest_weight = (--history_order.end());
            const uint64_t highest_id = highest_weight->second;
            const uint64_t highest_size = get<1>(history_map[highest_id]);
            if(current_size + highest_size <= cache_size) {
                current_size += highest_size;
                cache_set.insert(highest_id);
            }
            history_order.erase(highest_weight);
        }
        //        history_map.clear();
        cerr << "history activated - cache size: " << double(current_size)/cache_size << "\n";
    }
};




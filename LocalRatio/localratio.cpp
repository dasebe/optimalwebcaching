#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <stack>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>

// uncomment to enable debugging:
//#define CDEBUG 1
//#define EDEBUG 1

// util for debug
#ifdef CDEBUG
#define LOGnl(m,x,y,z,a) log_message(m,x,y,z,a,"\n")
#define LOG(m,x,y,z,a) log_message(m,x,y,z,a,"")
#else
#define LOGnl(m,x,y,z,a)
#define LOG(m,x,y,z,a)
#endif
#ifdef EDEBUG
#define LOGE(m,x,y,z,a) log_message(m,x,y,z,a,"\n")
#else
#define LOGE(m,x,y,z,a)
#endif
inline void log_message(std::string m, double x, double y, double z, double a, std::string e) {
    std::cerr << m << "," << x << "," << y  << "," << z << "," << a << e;
}

struct trEntry {
    const uint64_t id;
    const uint64_t size;
    double pval;
    size_t nextSeen;
    bool inSchedule;
    bool hasNext;

    trEntry(uint64_t nid, uint64_t nsize)
        : id(nid),
          size(nsize)
    {
        pval = 1;
        nextSeen = 0;
        inSchedule = false;
        hasNext = false;
    };
};

int main(int argc, char* argv[]) {

    // parameters
    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));

    // each interval is represented by the index of the trace, where it started
    // store trace here
    std::vector<trEntry> trace;
    // store lastSeen entries here
    std::map<std::pair<uint64_t, uint64_t>, size_t> lastSeen;
    // store which intervals are in the schedule
    std::unordered_set< size_t > scheduleSet;

    // open trace file
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0;
    // scan trace and initialize data structures
    while(traceFile >> time >> id >> size) {
        // only consider objects that are requested at least once (ignore last request, because doesn't make sense to cache that one)
        if(lastSeen.count(std::make_pair(id,size))>0) {
            // see object second time: add to schedule and initialize datastructures
            // find trace index where this interval started
            const size_t indexLastSeen = lastSeen[std::make_pair(id,size)];
            // add to schedule the start index of this interval
            // TODO MIGHT NOT NEED THIS
            scheduleSet.insert(indexLastSeen);
            trEntry & lastSeenEntry = trace[indexLastSeen];
            // update trace so that we know this interval has finite length
            lastSeenEntry.hasNext = true;
            // update trace so that we know the end time of this interval
            lastSeenEntry.nextSeen = reqc;
            // update trace: this interval is in schedule
            lastSeenEntry.inSchedule = true;
        }
        // store: id, size, bool seen before, pvalue, start of this interval (if started)
        trace.emplace_back(id,size);
        lastSeen[std::make_pair(id,size)]=reqc++;
    }

    // free memory
    lastSeen.clear();
    std::cerr << "parsed trace\n";
    
    // DEBUG print trace
#ifdef CDEBUG
    reqc = 0;
    LOGnl("TRACE",0,0,0,0);
    for(auto & curEntry : trace) {
        LOGnl(" ",curEntry.id,curEntry.size,0,curEntry.nextSeen);
    }
#endif
    // END DEBUG
    

    // create feasible schedule by excluding worst feasibility violations 
    // i.e., select largest delta in every iteration

    // store excluded intervals in a stack for later
    std::stack<size_t> excluded;

    // delta data structures
    int64_t curDelta;
    int64_t deltaStar;
    size_t timeStar;
    // map: intersecting interval set
    std::map<std::pair<uint64_t, uint64_t>, size_t> curI; //intersecting intervals at current time
    std::map<std::pair<uint64_t, uint64_t>, size_t> maxI; // intersectin intervals at timeStar

    // iterate (every iteration removes one interval, so max iteration count = trace size)
    for(size_t i=0; i<trace.size(); i++) {

        // find highest delta
        reqc = 0;
        // iterate over all time instances (we can skip the last one)

        maxI.clear();
        curDelta = -cacheSize;
        deltaStar = -cacheSize;
        timeStar = 0;
        for(size_t j=0; j<trace.size(); j++) {
            trEntry & curEntry = trace[j];

            // if no next request and in intersecting intervals -> remove
            if(!curEntry.hasNext & curI.count(std::make_pair(curEntry.id,curEntry.size))>0) {
                curI.erase(std::make_pair(curEntry.id,curEntry.size));
                curDelta -= curEntry.size;
                assert(!curEntry.inSchedule);
            }

            // if in schedule -> update curI
            if(curEntry.inSchedule) {

                // if not already in current intersecting set
                if(curI.count(std::make_pair(curEntry.id,curEntry.size))<=0 ) {
                    curDelta += curEntry.size;
                } // else: don't need update the size/width

                // add to current intersecting set
                curI[std::make_pair(curEntry.id,curEntry.size)] = j;

                // check if we need to update deltaStar
                if(curDelta > deltaStar) {
                    deltaStar = curDelta;
                    timeStar = j;
                    // copy intersecting set
                    maxI = curI;
                }
            }

#ifdef CDEBUG
            for(auto it: curI) {
                trEntry & curEntry = trace[it.second];
                LOG("|",curEntry.id,curEntry.size,0,0);
            }
            LOGnl("|  TW ",curDelta,j,deltaStar,timeStar);
#endif                
        }

        // check if we found a feasible solution
        if(deltaStar <= 0)
            break;

        curI.clear();
        assert(maxI.size()>0);
        assert(deltaStar > 0);
        assert(timeStar > 0); // first time interval can't be the unstable one
        
        LOGnl("\nd*,t*",deltaStar,timeStar,0,0);
        
        // find smallest rho so that p2 reaches zero
        double rho = 1;
        size_t toExclude = 0;
        for(auto & vit : maxI) {
            trEntry & curEntry = trace[vit.second];
            const double thisRho = curEntry.pval/static_cast<double>(curEntry.size);
            if(thisRho <= rho) {
                rho = thisRho;
                toExclude = vit.second;
                LOGE("mrho",rho,vit.second,curEntry.pval,curEntry.size);
            }
        }
        assert(rho < 1);
        LOGnl("min rho ",rho,0,0,0);
        
        // update p2, exclude intervals with p2=0 from schedule
        for(auto & vit : maxI) {
            trEntry & curEntry = trace[vit.second];

            // if (min-rho-entry or p2==0): we need the first check due to double rounding errors
            if(vit.second==toExclude || curEntry.pval <= rho * curEntry.size ) {
                LOGnl("exclude (t,id,size,p2)",vit.second,id,curEntry.size,curEntry.pval - rho*curEntry.size);
                curEntry.pval = 0;
                curEntry.inSchedule = false;
                // add current interval to excluded ones
                excluded.push(vit.second);
                // delete from current schedule
                scheduleSet.erase(vit.second);
            } else {
                // p2 >= 0
                curEntry.pval -= curEntry.size <= deltaStar ? rho * curEntry.size : rho * deltaStar; // pval = pval - p1
            }
        }
    }
    maxI.clear();
    std::cerr << "found feasible schedule\n";
    std::cout << "LR1 " << cacheSize << " hitc " << scheduleSet.size() << " reqc " << trace.size() << " OHR " << static_cast<double>(scheduleSet.size())/trace.size() << "\n";

    // we now have a feasible schedule, which might not be maximal

    while (!excluded.empty())
    {
        const size_t firstSeen = excluded.top();
        trEntry & newEntry = trace[excluded.top()];
        curI.clear();
        curDelta = -cacheSize;
        deltaStar = -cacheSize;

        for(size_t j=0; j<trace.size(); j++) {
            trEntry & curEntry = trace[j];

            // if no next request and in intersecting intervals -> remove
            if(!curEntry.hasNext & curI.count(std::make_pair(curEntry.id,curEntry.size))>0) {
                curI.erase(std::make_pair(curEntry.id,curEntry.size));
                curDelta -= curEntry.size;
                assert(!curEntry.inSchedule);
            }

            // if in schedule -> update curI
            if(j==excluded.top() || curEntry.inSchedule) {

                // if not already in current intersecting set
                if(curI.count(std::make_pair(curEntry.id,curEntry.size))<=0 ) {
                    curDelta += curEntry.size;
                } // else: don't need update the size/width

                // add to current intersecting set
                curI[std::make_pair(curEntry.id,curEntry.size)] = j;

                // check if we need to update deltaStar
                if(curDelta > deltaStar) {
                    deltaStar = curDelta;
                }
            }
        }
        // if still feasible add excluded.top() to schedule
        if(deltaStar <=0) {
            newEntry.inSchedule = true;
            scheduleSet.insert(excluded.top());
        }

        // pop stack
        excluded.pop();
    }

    std::cerr << "found maximal schedule\n";
    std::cout << "LR2 " << cacheSize << " hitc " << scheduleSet.size() << " reqc " << trace.size() << " OHR " << static_cast<double>(scheduleSet.size())/trace.size() << "\n";
    
    return 0;
}

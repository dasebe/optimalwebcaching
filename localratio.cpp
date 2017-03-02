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

// util for debug
#ifdef CDEBUG
#define LOGnl(m,x,y,z,a) log_message(m,x,y,z,a,"\n")
#define LOG(m,x,y,z,a) log_message(m,x,y,z,a,"")
#else
#define LOGnl(m,x,y,z,a)
#define LOG(m,x,y,z,a)
#endif
inline void log_message(std::string m, double x, double y, double z, double a, std::string e) {
    std::cerr << m << "," << x << "," << y  << "," << z << "," << a << e;
}


typedef std::tuple<uint64_t,uint64_t,bool,double,size_t> traceEntry; //id, size, hasNext, pvalue, nextSeen (index of next request)

int main(int argc, char* argv[]) {

    // parameters
    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));

    // each interval is represented by the index of the trace, where it started
    // store trace here
    std::vector<traceEntry> trace;
    // store lastSeen entries here
    std::map<std::pair<uint64_t, uint64_t>, size_t> lastSeen;
    // intervals that intersect at a time, point to request index
    std::vector<std::map< std::pair<uint64_t, uint64_t>, size_t > > sameTime;
    // store which intervals are in the schedule
    std::unordered_set< size_t > inSchedule;

    // open trace file
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    // scan trace and initialize data structures
    while(traceFile >> time >> id >> size) {
        // only consider objects that are requested at least once (ignore last request, because doesn't make sense to cache that one)
        if(lastSeen.count(std::make_pair(id,size))>0) {
            // see object second time: add to schedule and initialize datastructures
            // find trace index where this interval started
            const size_t indexLastSeen = lastSeen[std::make_pair(id,size)];
            // add to schedule the start index of this interval
            inSchedule.insert(indexLastSeen);
            // update trace so that we know this interval has finite length
            std::get<2>(trace[indexLastSeen]) = true;
            // update trace so that we know the end time of this interval
            std::get<4>(trace[indexLastSeen]) = reqc;
            // intersection set of this interval:
            // store indexLastSeen for every discrete time step this request has spanned
            for(uint64_t i=indexLastSeen; i<reqc; i++) {
                (sameTime[i])[std::make_pair(id,size)] = indexLastSeen;
            }
        } else {
            // see for first time, so count unique objects
            uniqc++;
        }
        // store: id, size, bool seen before, pvalue, start of this interval (if started)
        trace.emplace_back(id,size,false,1,0);
        sameTime.push_back( std::map< std::pair<uint64_t, uint64_t>, size_t >() );
        lastSeen[std::make_pair(id,size)]=reqc++;
    }

    std::cerr << "parsed trace\n";
    
    // DEBUG print trace
    reqc = 0;
    LOGnl("TRACE",0,0,0,0);
    for(traceEntry & thisTrEntry: trace) {
        const uint64_t id = std::get<0>(thisTrEntry);
        const uint64_t size = std::get<1>(thisTrEntry);
        double & pval = std::get<3>(thisTrEntry);
        const size_t firstSeen = reqc++;
        const size_t nextSeen = std::get<4>(thisTrEntry);
        LOGnl(" ",id,size,firstSeen,nextSeen);
    }
    // END DEBUG
    
    // initialize time width (total sum of object sizes that span a discrete time step)
    std::multimap<int64_t, uint64_t > widthToTime;
    std::unordered_map<uint64_t, std::multimap<int64_t, uint64_t >::iterator > timeToWidthIterator; // reverse mapping
    uint64_t currentWidth;
    reqc = 0;
    for(auto & it: sameTime) {
        currentWidth = 0;
        for(auto & vit: it) {
            const traceEntry & thisTrEntry = trace[vit.second];
            const uint64_t id = std::get<0>(thisTrEntry);
            const uint64_t size = std::get<1>(thisTrEntry);
            LOG("|",id,size,0,0);
            currentWidth+=size;
        }
        const int64_t currentDelta = currentWidth-cacheSize;
        auto tWIt = widthToTime.emplace(currentDelta, reqc);
        timeToWidthIterator[reqc] = tWIt;
        LOGnl("|  TW ",currentDelta,reqc,0,0);
        reqc++;
    }


    // create feasible schedule by excluding worst feasibility violations 
    // i.e., select largest delta in every iteration

    // store excluded intervals in a stack for later
    std::stack<size_t> excluded;

    // find largest delta
    auto widthIt = --widthToTime.end();
    int64_t deltaStar = widthIt->first;

    // iterate
    while (deltaStar > 0) {

        const uint64_t timeStar = widthIt->second;
        LOGnl("\nd*,t* ",deltaStar,timeStar,0,0);
        
        // find smallest rho so that p2 reaches zero
        double rho = 1;
        for(auto & it : sameTime[timeStar]) {
            const traceEntry & thisTrEntry = trace[it.second];
            const uint64_t id = std::get<0>(thisTrEntry);
            const double size = std::get<1>(thisTrEntry);
            const double pval = std::get<3>(thisTrEntry);
            const double thisRho = pval/size;
            rho = thisRho <= rho ? thisRho : rho;
        }
        assert(rho < 1);
        LOGnl("min rho ",rho,0,0,0);

        // update p2, exclude p2=0 intervals from schedule
        for(auto it=sameTime[timeStar].begin(); it!=sameTime[timeStar].end(); ) {

            // define constants with this object's properties
            traceEntry & thisTrEntry = trace[it->second];
            const uint64_t id = std::get<0>(thisTrEntry);
            const uint64_t size = std::get<1>(thisTrEntry);
            double & pval = std::get<3>(thisTrEntry);
            const size_t firstSeen = it->second;
            const size_t nextSeen = std::get<4>(thisTrEntry);

            // p2==0
            if(pval <= rho * size ) {
                LOGnl("exclude (t,id,size,p2)",it->second,id,size,pval - rho*size);

                // exclude from schedule and update width  of all times between firstSeen and nextSeen
                for(size_t i=firstSeen; i<nextSeen; i++) {
                    // find iterator to widthToTime map and assert it exists
                    auto wTITold = timeToWidthIterator.find(i);
                    assert(wTITold != timeToWidthIterator.end());
                    auto tWItold = timeToWidthIterator[i]; // iterator to this widthToTime field
                    assert(i==tWItold->second);
                    // calculate new delta
                    const int64_t oldDelta = tWItold->first;
                    const int64_t newDelta = oldDelta-size;
                    // delete from widthToTime map
                    widthToTime.erase(tWItold);
                    // create new entry in widthToTime with new delta value
                    auto tWItnew = widthToTime.emplace(newDelta, i);
                    timeToWidthIterator[i] = tWItnew;
                    // skip current it iterator
                    if(i!=timeStar) {
                        (sameTime[i]).erase(std::make_pair(id,size));
                    }
                    // added to queue of excluded intervals
                }
                // set p2 to 0 by definition
                pval = 0;
                // add current interval to excluded ones
                excluded.push(it->second);
                // delete from current schedule
                inSchedule.erase(it->second);
                // delete current entry from current sameTime map
                it = sameTime[timeStar].erase(it);
            } else {
                // p2 >= 0
                pval -= size <= deltaStar ? rho * size : rho * deltaStar; // pval = pval - p1
                it++; // continue iterator
            }
        }
        // find time step with max delta
        assert(widthToTime.size()>0);

        //debug
        reqc = 0;
        for(auto & it: sameTime) {
            for(auto & vit: it) {
                const traceEntry & thisTrEntry = trace[vit.second];
                const uint64_t id = std::get<0>(thisTrEntry);
                const uint64_t size = std::get<1>(thisTrEntry);
                LOG("|",id,size,0,0);
            }
            LOGnl("|  TW ",timeToWidthIterator[reqc]->first,reqc,0,0);
            reqc++;
        }
        // end of debug
        
        widthIt = --widthToTime.end();
        deltaStar = widthIt->first;
    }

    std::cerr << "found feasible schedule\n";

    // we now have a feasible schedule, which might not be maximal

    while (!excluded.empty())
    {
        traceEntry & thisTrEntry = trace[excluded.top()];
        const uint64_t id = std::get<0>(thisTrEntry);
        const uint64_t size = std::get<1>(thisTrEntry);
        const size_t firstSeen = excluded.top();
        const size_t nextSeen = std::get<4>(thisTrEntry);
        bool feasibleToAdd = true;

        // check consistency (shouldn't be necessary)
        assert(inSchedule.find(firstSeen)==inSchedule.end());

        // check feasibility for all time steps
        for(size_t i=firstSeen; i<nextSeen; i++) {
            // find iterator to widthToTime map and assert it exists
            auto wTITold = timeToWidthIterator.find(i);
            assert(wTITold != timeToWidthIterator.end());
            auto tWItold = timeToWidthIterator[i]; // iterator to this widthToTime field
            assert(i==tWItold->second);
            // calculate new delta
            const int64_t oldDelta = tWItold->first;
            const int64_t newDelta = oldDelta+size;
            if(newDelta > 0) {
                feasibleToAdd = false;
                break;
            }
        }

        // actually add this one
        if(feasibleToAdd) {
            LOGnl("adding to schedule",firstSeen,0,0,0);
            // add to feasible schedule
            inSchedule.insert(firstSeen);
            // update intersecting time steps
            for(size_t i=firstSeen; i<nextSeen; i++) {
                // find iterator to widthToTime map and assert it exists
                auto wTITold = timeToWidthIterator.find(i);
                assert(wTITold != timeToWidthIterator.end());
                auto tWItold = timeToWidthIterator[i]; // iterator to this widthToTime field
                assert(i==tWItold->second);
                // calculate new delta
                const int64_t oldDelta = tWItold->first;
                const int64_t newDelta = oldDelta+size;
                // delete from widthToTime map
                widthToTime.erase(tWItold);
                // create new entry in widthToTime with new delta value
                auto tWItnew = widthToTime.emplace(newDelta, i);
                timeToWidthIterator[i] = tWItnew;
            }
        }
        
        excluded.pop();
    }

    std::cerr << "found maximal schedule\n";

    std::cout << "hitc " << inSchedule.size() << " reqc " << trace.size() << " OHR " << static_cast<double>(inSchedule.size())/trace.size() << "\n";
    
    return 0;
}

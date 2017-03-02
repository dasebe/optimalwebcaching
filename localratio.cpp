#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <functional>

typedef std::tuple<uint64_t,uint64_t,bool,double,size_t> traceEntry; //id, size, hasNext, pvalue, nextSeen (index of next request)

int main(int argc, char* argv[]) {

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    
    // parse trace file
    std::vector<traceEntry> trace;
    std::map<std::pair<uint64_t, uint64_t>, size_t> lastSeen;
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    std::vector<std::map< std::pair<uint64_t, uint64_t>, size_t > > sameTime; // intervals that intersect at a time, point to request index

    // scan trace and initialize data structures
    while(traceFile >> time >> id >> size) {
        // only consider objects that are requested at least once (ignore last request, because doesn't make sense to cache that one)
        if(lastSeen.count(std::make_pair(id,size))>0) {
            // see object second time
            const size_t indexLastSeen = lastSeen[std::make_pair(id,size)];
            std::get<2>(trace[indexLastSeen]) = true;
            std::get<4>(trace[indexLastSeen]) = reqc;
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

    // DEBUG print trace
    reqc = 0;
    std::cout << " TRACE \n";
    for(traceEntry & thisTrEntry: trace) {
        const uint64_t id = std::get<0>(thisTrEntry);
        const uint64_t size = std::get<1>(thisTrEntry);
        double & pval = std::get<3>(thisTrEntry);
        const size_t firstSeen = reqc++;
        const size_t nextSeen = std::get<4>(thisTrEntry);
        std::cout << id << " " << size << " " << pval << " " << firstSeen << " " << nextSeen << "\n";
    }
    std::cout << "\n\n";
    // END DEBUG
    
    // initialize time width (total sum of object sizes that span a discrete time step)
    std::multimap<uint64_t, uint64_t > widthToTime;
    std::unordered_map<uint64_t, std::multimap<uint64_t, uint64_t >::iterator > timeToWidthIterator; // reverse mapping
    uint64_t thisWidthToTime;
    reqc = 0;
    for(auto & it: sameTime) {
        thisWidthToTime = 0;
        for(auto & vit: it) {
            const traceEntry & thisTrEntry = trace[vit.second];
            const uint64_t id = std::get<0>(thisTrEntry);
            const uint64_t size = std::get<1>(thisTrEntry);
            std::cout << "(" << id << "," << size << ") ";
            thisWidthToTime+=size;
        }
        auto tWIt = widthToTime.emplace(thisWidthToTime, reqc);
        timeToWidthIterator[reqc] = tWIt;
        std::cout << " TW " << thisWidthToTime << " " << reqc;
        std::cout << "\n";
        reqc++;
    }

    // largest width
    auto widthIt = --widthToTime.end();
    uint64_t deltaStar = widthIt->first;


    while (deltaStar > 0) {

        const uint64_t timeStar = widthIt->second;
        std::cout <<  "\n\n\nd* " << deltaStar << " t* " <<  timeStar<< "\n";
        // find request interval with largest size
        uint64_t maxSize = 0;
        for(auto & it : sameTime[timeStar]) {
            const traceEntry & thisTrEntry = trace[it.second];
            const uint64_t id = std::get<0>(thisTrEntry);
            const uint64_t size = std::get<1>(thisTrEntry);
            maxSize = size > maxSize ? size : maxSize;
            std::cout << "(" << id << "," << size << ") | ";
        }
        std::cout << "\n";
        assert(maxSize > 0);
        const double rho = maxSize <= deltaStar ? 1.0/maxSize : 1.0/deltaStar;
        std::cout << "maxSize " << maxSize << " rho " << rho << "\n";
        for(auto it=sameTime[timeStar].begin(); it!=sameTime[timeStar].end(); ) {

            // define constants with this object's properties
            traceEntry & thisTrEntry = trace[it->second];
            const uint64_t id = std::get<0>(thisTrEntry);
            const uint64_t size = std::get<1>(thisTrEntry);
            double & pval = std::get<3>(thisTrEntry);
            const size_t firstSeen = it->second;
            const size_t nextSeen = std::get<4>(thisTrEntry);

            // check if this p2 function is zero
            if(pval <= rho * size ) {
                std::cout << "throw out (" << id << "," << size << ") pval " << pval << " dS " << rho*size << " p2 " << pval - rho*size << "\n";
                // update width of all times between firstSeen and nextSeen
                // store indexLastSeen for every discrete time step this request has spanned
                for(uint64_t i=firstSeen; i<nextSeen; i++) {
                    auto wTITold = timeToWidthIterator.find(i);
                    assert(wTITold != timeToWidthIterator.end());
                    auto tWItold = timeToWidthIterator[i]; // iterator to this widthToTime field
                    assert(i==tWItold->second);
                    const uint64_t delta = tWItold->first;
                    widthToTime.erase(tWItold);
                    auto tWItnew = widthToTime.emplace(delta-size, i);
                    timeToWidthIterator[i] = tWItnew;
                    // skip current it iterator
                    if(i!=timeStar) {
                        (sameTime[i]).erase(std::make_pair(id,size));
                    }
                }
                pval = 0; // p2 = zero
                it = sameTime[timeStar].erase(it); // delete it and continue iterator
            } else {
                pval -= size <= deltaStar ? rho * size : rho * deltaStar; // p2>= zero
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
                std::cout << "(" << id << "," << size << ") ";
            }
            std::cout << " TW " << timeToWidthIterator[reqc]->first << " " << reqc;
            std::cout << "\n";
            reqc++;
        }
        // end of debug
        
        widthIt = --widthToTime.end();
        deltaStar = widthIt->first;
    }
    return 0;
}

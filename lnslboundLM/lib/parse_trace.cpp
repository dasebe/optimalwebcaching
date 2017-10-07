//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <tuple>
#include <cmath>
#include <unordered_set>
#include "parse_trace.h"

//using namespace lemon;

uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path,uint64_t cacheSize) {
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;
    std::unordered_set<std::pair<uint64_t, uint64_t> > cacheable;

    while(traceFile >> time >> id >> size) {
        const auto idSize = std::make_pair(id,size);
        if(size <= cacheSize) {
            lastSeen[idSize]++;
        }
        reqc++;
    }
       
    for(auto it: lastSeen) {
        if(it.second>1) {
            cacheable.insert(it.first);
        }
    }

    //reset
    std::cerr << "lines: " << reqc << " objects " << lastSeen.size() << "\n";
    trace.reserve(reqc - lastSeen.size());
    lastSeen.clear();
    reqc=0;
    traceFile.clear();
    traceFile.seekg(0, std::ios::beg);

    while(traceFile >> time >> id >> size) {
        const auto idSize = std::make_pair(id,size);
        if(lastSeen.count(idSize)>0) {
            trace[lastSeen[idSize]].nextSeen = reqc;
            const double intervalLength = reqc-lastSeen[idSize];
            // calculate utility
            const double utilityDenominator = size*intervalLength;
            assert(utilityDenominator>0);
            trace[lastSeen[idSize]].utility = 1.0L/utilityDenominator;
        } else {
            uniqc++;
        }
        // only add if object is cacheable
        if(cacheable.count(idSize)>0) {
            trace.emplace_back(id,size);
            lastSeen[idSize]=reqc++;
        }
    }
    return uniqc;
}
                    
    uint64_t createMCF(SmartDigraph & g, std::vector<trEntry > & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int64_t> & supplies, const double minUtil, const double maxUtil, std::unordered_map<size_t, int> & traceToArc) {

    size_t effectiveEjectSize = 0;

    // create a graph with just arcs with utility between minUtil and maxUtil
    // mcf instance data
    SmartDigraph::Node curNode = g.addNode(); // initial node
    SmartDigraph::Arc curArc;
    SmartDigraph::Node prevNode;

    // track cached/non-cached intervals
    std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, int> > lastSeen;
    long double nonFlexSize = 0;
    std::map<size_t, long double> endOfIntervalSize;


    for(uint64_t i=0; i<trace.size(); i++) {
        trEntry & curEntry = trace[i];

        // first: check if previous interval ended here
        if(lastSeen.count(std::make_pair(curEntry.id,curEntry.size))>0) {
            // create "outer" request arc
            const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[std::make_pair(curEntry.id,curEntry.size)].second);
            curArc = g.addArc(lastReq,curNode);
            cap[curArc] = curEntry.size;
            cost[curArc] = 1/static_cast <double>(curEntry.size);
            supplies[lastReq] += curEntry.size;
            supplies[curNode] -= curEntry.size;
            traceToArc[lastSeen[std::make_pair(curEntry.id,curEntry.size)].first] = g.id(curArc);
            effectiveEjectSize++;
            // DEBUG
            LOG("oA",i,lastSeen[std::make_pair(curEntry.id,curEntry.size)].first,g.id(curArc));
            // delete lastSeen entry, as the next interval might not be part
            lastSeen.erase(std::make_pair(curEntry.id,curEntry.size));
        } 

        // create arcs if in ejection set
        if(isInEjectSet(minUtil, maxUtil, curEntry) ) {
            // second: if there is another request for this object
            if(curEntry.nextSeen > 0) {
                // save prev node as anchor for future arcs
                prevNode = curNode;
                lastSeen[std::make_pair(curEntry.id,curEntry.size)]=std::make_pair(i,g.id(prevNode));
                // create another node, "inner" capacity arc
                curNode = g.addNode(); // next node
                curArc = g.addArc(prevNode,curNode);
                cap[curArc] = cacheSize - std::floor(nonFlexSize); 
                // DEBUG
                LOG("iA",i,cap[curArc],0);
                cost[curArc] = 0;
            }

            //not in ejection set and dvar > 0 -> need to subtract dvar*size for interval's duration
        } else if (curEntry.dvar > 0) { 
            const long double curEffectiveSize = curEntry.size*curEntry.dvar;
            // DEBUG
            LOG("-CS",i,curEffectiveSize,curEntry.dvar);
            assert(cacheSize >= curEffectiveSize);
            nonFlexSize += curEffectiveSize;
            // assert valid flexsize
            endOfIntervalSize.emplace(curEntry.nextSeen,curEffectiveSize);
        }


        // clear all nonFlexSize which are 
        while(endOfIntervalSize.size() > 0 && endOfIntervalSize.begin()->first <= i+1) {
            // DEBUG
            LOG("+CS",i,endOfIntervalSize.begin()->first,endOfIntervalSize.begin()->second);
            nonFlexSize -= endOfIntervalSize.begin()->second;
            endOfIntervalSize.erase(endOfIntervalSize.begin());
        }
    }

    return effectiveEjectSize;

}



bool feasibleCacheAll(std::vector<trEntry > & trace, uint64_t cacheSize, const double minUtil) {

    // delta data structures: from localratio technique
    int64_t curDelta;
    int64_t deltaStar;
    // map currently cached (id,size) to interval index in trace
    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t> curI; //intersecting intervals at current time

    curDelta = -cacheSize;
    deltaStar = -cacheSize;

    for(size_t j=0; j<trace.size(); j++) {
        trEntry & curEntry = trace[j];

        // if no next request and in intersecting intervals -> remove
        if( (curEntry.nextSeen == 0) & (curI.count(std::make_pair(curEntry.id,curEntry.size)) > 0) ) {
            curI.erase(std::make_pair(curEntry.id,curEntry.size));
            curDelta -= curEntry.size;
            assert(curEntry.dvar==0);
        }

        // if with utility in [minUtil,1]
        if(isInEjectSet(minUtil, 1.01, curEntry) && cacheSize >= curEntry.size) {

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
    // return feasibility bool
    return (deltaStar <=0);

}

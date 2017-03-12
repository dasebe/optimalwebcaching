//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include "parse_trace.h"

//using namespace lemon;

uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path) {
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;

    while(traceFile >> time >> id >> size) {
        if(lastSeen.count(std::make_pair(id,size))>0) {
            trace[lastSeen[std::make_pair(id,size)]].hasNext = true;
            trace[lastSeen[std::make_pair(id,size)]].nextSeen = reqc;
            const long double intervalLength = reqc-lastSeen[std::make_pair(id,size)];
            // calculate utility
            const long double utilityDenominator = size*intervalLength;
            assert(utilityDenominator>0);
            trace[lastSeen[std::make_pair(id,size)]].utility = 1.0L/utilityDenominator;
        } else {
            uniqc++;
        }
        trace.emplace_back(id,size,time);
        lastSeen[std::make_pair(id,size)]=reqc++;
    }
    return uniqc;
}
                    
uint64_t createMCF(SmartDigraph & g, std::vector<trEntry > & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int64_t> & supplies, const long double minUtil, const long double maxUtil) {

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
        curEntry.active = false;

        // first: check if previous interval ended here
        if(lastSeen.count(std::make_pair(curEntry.id,curEntry.size))>0) {
            // create "outer" request arc
            const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[std::make_pair(curEntry.id,curEntry.size)].second);
            curArc = g.addArc(lastReq,curNode);
            cap[curArc] = curEntry.size;
            cost[curArc] = 1/static_cast <double>(curEntry.size);
            supplies[lastReq] += curEntry.size;
            supplies[curNode] -= curEntry.size;
            trace[lastSeen[std::make_pair(curEntry.id,curEntry.size)].first].arcId = g.id(curArc);
            trace[lastSeen[std::make_pair(curEntry.id,curEntry.size)].first].active = true;
            effectiveEjectSize++;
            // DEBUG
            LOG("oA",i,lastSeen[std::make_pair(curEntry.id,curEntry.size)].first,g.id(curArc));
            // delete lastSeen entry, as the next interval might not be part
            lastSeen.erase(std::make_pair(curEntry.id,curEntry.size));
        } 

        // create arcs if in ejection set
        if(isInEjectSet(minUtil, maxUtil, curEntry) ) {
            // second: if there is another request for this object
            if(curEntry.hasNext) {
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

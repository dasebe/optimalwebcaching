#include <fstream>
#include <cassert>
#include <map>
#include <unordered_map>
#include <tuple>
#include "parse_trace.h"

using namespace lemon;

uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path) {
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    std::map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;

    while(traceFile >> time >> id >> size) {
        if(lastSeen.count(std::make_pair(id,size))>0) {
            trace[lastSeen[std::make_pair(id,size)]].hasNext = true;
            const long double intervalLength = reqc-lastSeen[std::make_pair(id,size)];
            assert(intervalLength>0);
            trace[lastSeen[std::make_pair(id,size)]].iLen = intervalLength;
        } else {
            uniqc++;
        }
        trace.emplace_back(id,size,time);
        lastSeen[std::make_pair(id,size)]=reqc++;
    }
    return uniqc;
}

SmartDigraph::Node createMCF(SmartDigraph & g, std::vector<trEntry> & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<long double> & cost, SmartDigraph::NodeMap<int64_t> & supplies) {

    // we consider (id,size) as unique identification of an object (sizes can change, but then it's a different object)
    // lastSeen maps (id,size) to (nodeId,traceIndex) of the last time this object was seen
    std::unordered_map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, int> > lastSeen;
    SmartDigraph::Arc curInnerArc, curOuterArc;
    SmartDigraph::Node curNode = g.addNode(); // initial node
    SmartDigraph::Node prevNode;

    // iterate over trace
    for(uint64_t i=0; i<trace.size(); i++) {
        trEntry & curEntry = trace[i];
        const uint64_t id = curEntry.id;
        const uint64_t size = curEntry.size;
        // first: check if previous interval ended here
        if(lastSeen.count(std::make_pair(id,size))>0) {
            // create "outer" request arc
            const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[std::make_pair(id,size)].second);
            curOuterArc = g.addArc(lastReq,curNode);
            cap[curOuterArc] = size;
            cost[curOuterArc] = 1/static_cast <long double>(size);
            supplies[lastReq] += size;
            supplies[curNode] -= size;
            lastSeen.erase(std::make_pair(id,size));
        }
        // second: if there is another request for this object
        if(curEntry.hasNext) {
            // save prev node as anchor for future arcs
            prevNode = curNode;
            lastSeen[std::make_pair(id,size)]=std::make_pair(i,g.id(prevNode));
            // create another node, "inner" capacity arc
            curNode = g.addNode(); // next node
            curInnerArc = g.addArc(prevNode,curNode);
            cap[curInnerArc] = cacheSize; 
            cost[curInnerArc] = 0;
        }
        // set current inner arc id
        curEntry.innerArcId = g.id(curInnerArc);
    }

    assert(lastSeen.size()==0);

    return curNode;
}

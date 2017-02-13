#include <iostream>
#include <cassert>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>
#include <lemon/smart_graph.h>

using namespace lemon;

void parseTraceFile(std::vector<std::tuple<uint64_t,uint64_t,bool> > & trace, std::string & path) {
    std::ifstream traceFile(path);
    uint64_t t, id, size, reqc=0, uniqc=0;
    std::map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;

    while(traceFile >> t >> id >> size) {
        if(lastSeen.count(std::make_pair(id,size))>0) {
            std::get<2>(trace[lastSeen[std::make_pair(id,size)]]) = true;
        } else {
            uniqc++;
        }
        trace.emplace_back(id,size,false);
        lastSeen[std::make_pair(id,size)]=reqc++;
    }
    std::cout << "scanned trace n=" << reqc << " m= " << uniqc << std::endl;
}
                    
void createMCF(SmartDigraph & g, std::vector<std::tuple<uint64_t,uint64_t,bool> > & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int> & supplies) {

    // init
    std::map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;
    uint64_t reqc=0;
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curNode = g.addNode(); // initial node
    SmartDigraph::Node prevNode;

    // iterate over trace
    for(auto it: trace) {
        // only consider requests that reoccur
        // should we delete if does not reoccur?
        if(std::get<2>(it)) {
            const uint64_t id=std::get<0>(it);
            const uint64_t size=std::get<1>(it);
            if(lastSeen.count(std::make_pair(id,size))>0) {
                // create request arc ("outer")
                const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[std::make_pair(id,size)]);
                curArc = g.addArc(lastReq,curNode);
                cap[curArc] = size;
                cost[curArc] = 1/static_cast <double>(size);
                supplies[lastReq] += size;
                supplies[curNode] -= size;
            }
            // create capacity arc ("inner")
            prevNode = curNode;
            curNode = g.addNode(); // next node
            curArc = g.addArc(prevNode,curNode);
            cap[curArc] = cacheSize; 
            cost[curArc] = 0;
            lastSeen[std::make_pair(id,size)]=reqc;
            reqc++;
        }
    }
    
    // there will m intervals without arcs
    // they all intersect and thus point to the same final node
    for(auto it: lastSeen) {
        const uint64_t id=it.first.first;
        const uint64_t reqc=it.second;
        const uint64_t size=it.first.second;
        const SmartDigraph::Node lastReq = g.nodeFromId(reqc);
        curArc = g.addArc(lastReq,curNode);
        cap[curArc] = size;
        cost[curArc] = 1/static_cast <double>(size);
        supplies[lastReq] += size;
        supplies[curNode] -= size;
    }
}    

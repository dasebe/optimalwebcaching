//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <unordered_map>
#include <tuple>
#include "parse_trace.h"

using namespace lemon;

uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path) {
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;

    while(traceFile >> time >> id >> size) {
        if(lastSeen.count(std::make_pair(id,size))>0) {
            trace[lastSeen[std::make_pair(id,size)]].hasNext = true;
            trace[lastSeen[std::make_pair(id,size)]].nextSeen = reqc;
            trace[lastSeen[std::make_pair(id,size)]].dvar = 0;
        } else {
            uniqc++;
        }
        trace.emplace_back(id,size,time);
        lastSeen[std::make_pair(id,size)]=reqc++;
    }
    return uniqc;
}
                    
//void createMCF(SmartDigraph & g, std::vector<trEntry > & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int64_t> & supplies) {

    // we consider (id,size) as unique identification of an object (sizes can change, but then it's a different object)
    // lastSeen maps (id,size) to (nodeId,traceIndex) of the last time this object was seen
//}

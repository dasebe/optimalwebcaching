//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <tuple>
#include <cmath>
#include <unordered_set>
#include "parse_trace.h"

//using namespace lemon;

uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path, uint64_t & uniqCount) {
    std::ifstream traceFile(path);
    uint64_t time, id, reqc=0, totalreqc=0;
    int64_t size;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;
    std::unordered_set<std::pair<uint64_t, uint64_t> > cacheable;

    while(traceFile >> time >> id >> size) {
        const auto idSize = std::make_pair(id,size);
        if(size > 0) {
            lastSeen[idSize]++;
        }
        totalreqc++;
    }
       
    for(auto it: lastSeen) {
        if(it.second>1) {
            cacheable.insert(it.first);
        }
    }

    //reset
    std::cerr << "lines: " << totalreqc << " objects " << lastSeen.size() << "\n";
    trace.reserve(totalreqc - lastSeen.size());
    lastSeen.clear();
    reqc=0;
    traceFile.clear();
    traceFile.seekg(0, std::ios::beg);
    uniqCount = 0;

    while(traceFile >> time >> id >> size) {
        const auto idSize = std::make_pair(id,size);
        if(lastSeen.count(idSize)>0) {
            trace[lastSeen[idSize]].nextSeen = reqc;
        } else {
            if(cacheable.count(idSize)>0) {
                uniqCount++;
            }
        }
        // only add if object is cacheable
        if(cacheable.count(idSize)>0) {
            assert(size>0);
            trace.emplace_back(size);
            lastSeen[idSize]=reqc++;
        }
    }
    return totalreqc;
}

//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include "parse_trace.h"

void parseTraceFile(std::vector<trEntry> & trace, std::string & path) {
    std::ifstream traceFile(path);

    uint64_t time, id, size, reqc=0;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;

    while(traceFile >> time >> id >> size) {
        const auto idsize = std::make_pair(id,size);
        if(lastSeen.count(idsize)>0) {
            trace[lastSeen[idsize]].hasNext = true;
            const uint64_t volume = (reqc-lastSeen[idsize]) * size;
            trace[lastSeen[idsize]].volume = volume;
        }
        trace.emplace_back(id,size,time);
        lastSeen[idsize]=reqc++;
    }
}

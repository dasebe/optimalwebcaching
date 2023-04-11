//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <cmath>
#include "parse_trace.h"

void parseTraceFile(std::vector<trEntry> & trace, std::string & path, uint64_t warmup) {
    std::ifstream traceFile(path);

    uint64_t time, id, size, reqc=0;
    uint64_t count = 0;
    uint64_t count2 = 0;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;
    std::unordered_set<std::pair<uint64_t, uint64_t>> warmup_objs;

    while(traceFile >> time >> id >> size) {
	count++;
	const auto idsize = std::make_pair(id,size);
	if (count < warmup && size > 0) {
	     warmup_objs.insert(idsize);
	     continue;
	}
	count2++;
	trace.emplace_back(size);
	if(size > 0 && lastSeen.count(idsize)==0 && warmup_objs.count(idsize)>0) {
	    trace[reqc].hasNext = true;
	    const uint64_t volume = (reqc-0) * size;
	    trace[reqc].volume = volume;
	} else if(size > 0 && lastSeen.count(idsize)>0) {
	    trace[lastSeen[idsize]].hasNext = true;
            const uint64_t volume = (reqc-lastSeen[idsize]) * size;
            trace[lastSeen[idsize]].volume = volume;
        }
        lastSeen[idsize]=reqc++;
    }
    std::cerr << "total requests: " << count << " after warmup: " << count2 << "\n";
}

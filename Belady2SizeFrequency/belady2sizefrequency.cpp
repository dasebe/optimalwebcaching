#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << argv[0] << " traceFile cacheSize" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    uint64_t cacheSize(std::stoull(argv[2]));

    // parse trace file
    std::vector<trEntry> trace;
    parseTraceFile(trace, path);
    uint64_t totalReqc = trace.size();
    std::cerr << "scanned trace n=" << totalReqc << " cs " << cacheSize << std::endl;

    // get nextSeen indices
    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t> lastSeen;
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> frequency;
    for (size_t i = trace.size(); i--> 0 ;) {
        trEntry & cur = trace[i];
        if(lastSeen.count(std::make_pair(cur.id,cur.size)) > 0) {
            cur.hasNext = true;
            cur.nextSeen = lastSeen[std::make_pair(cur.id, cur.size)];
            cur.frequency = frequency[std::make_pair(cur.id, cur.size)];
        }
        lastSeen[std::make_pair(cur.id, cur.size)] = i;
        frequency[std::make_pair(cur.id, cur.size)]++;
    }

    // actual caching algorithm
    cacheAlg(trace, cacheSize);

    // print results
    printRes(trace, "Belady2SizeForward "+std::to_string(cacheSize));

    // backward
    LOG("\n--------------------------\n\n",0,0,0);

    // reset
    lastSeen.clear();
    frequency.clear();
    // compute backwards
    for(size_t i=0; i<trace.size(); i++) {
        trEntry & cur = trace[i];
        // reset
        cur.hasNext = true;
        cur.nextSeen = 0;
        cur.hit = 0;
        cur.frequency = 1;
        if(lastSeen.count(std::make_pair(cur.id,cur.size)) > 0) {
            cur.nextSeen = lastSeen[std::make_pair(cur.id, cur.size)];
            cur.frequency = frequency[std::make_pair(cur.id, cur.size)];
        }
        lastSeen[std::make_pair(cur.id, cur.size)] = i;
        frequency[std::make_pair(cur.id, cur.size)]++;
    }

    // actual caching algorithm
    cacheAlg(trace, cacheSize);

    // print results
    printRes(trace, "Belady2SizeBackward "+std::to_string(cacheSize));
    
    return 0;
}



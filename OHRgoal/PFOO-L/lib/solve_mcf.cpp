#include "solve_mcf.h"
#include <fstream>
#include <algorithm>
#include <iomanip>

void cacheAlg(std::vector<trEntry> & trace) {
    std::sort(trace.begin(), trace.end());
}

void printRes(std::vector<trEntry> & trace, std::string algName, uint64_t cacheSizeMax, std::ofstream * resultfile) {
    long double csize = 0;
    long double nextCsizePrint = 1;
    long double lcacheSizeMax = cacheSizeMax;
    long double ltraceSize = trace.size();
    uint64_t hitc = 0;
    uint64_t reqcDiff = 0;
    *resultfile << std::fixed << std::setprecision(4);
    // iterate over sorted trace
    for(auto it: trace) {
        if(nextCsizePrint > lcacheSizeMax)
            break;
        if(csize>=nextCsizePrint) {
            *resultfile << algName << " " << nextCsizePrint << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << " " << csize << " " << reqcDiff << "\n";
            nextCsizePrint*=2;
            reqcDiff=0;
        }
        // add up hits and used-up "fluid" cache size
        if(it.hasNext) {
            hitc++;
            csize += it.volume/ltraceSize;
            reqcDiff++;
        }
    }
    // fill in the gaps (if the cache fits the whole trace)
    // keep outputting the last hit count (as trace is too short)
    while(nextCsizePrint <= lcacheSizeMax) {
        std::cerr << "filling in gaps, trace too short\n";
        *resultfile << algName << " " << nextCsizePrint << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << " " << csize << " " << reqcDiff << "\n";
        nextCsizePrint*=2;
    }
            
    resultfile->close();
}

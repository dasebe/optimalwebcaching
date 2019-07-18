#include "solve_mcf.h"
#include <fstream>
#include <algorithm>

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
    for(auto it: trace) {
        if(csize> lcacheSizeMax)
            break;
        if(csize>=nextCsizePrint) {
            *resultfile << algName << " " << nextCsizePrint << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << " " << csize << " " << reqcDiff << "\n";
            //            std::cout << algName << " " << csize << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << "\n";
            nextCsizePrint*=2;
            reqcDiff=0;
        }
        if(it.hasNext) {
            hitc++;
            csize += it.volume/ltraceSize;
            reqcDiff++;
        }
    }
    while(nextCsizePrint < lcacheSizeMax/2) {
        *resultfile << algName << " " << nextCsizePrint << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << " " << csize << " " << reqcDiff << "\n";
        nextCsizePrint*=2;
    }
            
    resultfile->close();
}

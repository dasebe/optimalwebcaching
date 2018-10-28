#include "solve_mcf.h"
#include <iostream>
#include <fstream>
#include <algorithm>

void cacheAlg(std::vector<trEntry> & trace) {
    std::sort(trace.begin(), trace.end());
}

void printRes(std::vector<trEntry> & trace, uint64_t & byteSum, uint64_t cacheSize) {
    uint64_t totalCacheVolume = cacheSize * trace.size();
    uint64_t currentVolume = 0;
    uint64_t hitc = 0;
    uint64_t bytehitc = 0;
    for(auto it: trace) {
        if(currentVolume > totalCacheVolume)
            break;
        if(it.hasNext) {
            hitc++;
            bytehitc += it.size;
            currentVolume += it.volume;
        }
    }
    std::cout << "PFOO-L ohr " << double(hitc)/trace.size() << " bhr " << double(bytehitc)/byteSum << "\n";
}

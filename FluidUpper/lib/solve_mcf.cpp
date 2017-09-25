#include "solve_mcf.h"
#include <fstream>
#include <algorithm>

void cacheAlg(std::vector<trEntry> & trace) {
    std::sort(trace.begin(), trace.end());
}

void printRes(std::vector<trEntry> & trace, std::string algName, uint64_t cacheSizeMax, std::ofstream * resultfile) {
    uint64_t csize = 0;
    uint64_t nextCsizePrint = 1;
    uint64_t hitc = 0;
    for(auto it: trace) {
        csize += it.volume/trace.size();
        if(csize>cacheSizeMax)
            break;
        if(it.hasNext)
            hitc++;
        if(csize>nextCsizePrint) {
            *resultfile << algName << " " << csize << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << "\n";
            //            std::cout << algName << " " << csize << " " << hitc << " " << trace.size() << " " << (double)hitc/trace.size() << "\n";
            nextCsizePrint*=2;
        }
    }
    resultfile->close();
}

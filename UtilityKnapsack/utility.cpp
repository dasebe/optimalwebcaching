#include <fstream>
#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "lib/parse_trace.h"
#include <algorithm>
#include <set>

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << argv[0] << " traceFile" << std::endl;
        return 1;
    }

    std::string path(argv[1]);

    std::ifstream traceFile(path);
    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> counts;

    uint64_t reqc=0;
    uint64_t time, id, size;
    while(traceFile >> time >> id >> size) {
        const auto idSize=std::make_pair(id,size);
        counts[idSize]++;
        reqc++;
    }

    std::vector<trEntry> prios;
    for(auto it: counts) {
        prios.emplace_back(std::get<1>(it.first),
                           it.second/double(std::get<1>(it.first)),
                           it.second);
    }
    std::sort(prios.begin(), prios.end());
    uint64_t cacheSize = 1;
    uint64_t curSize = 0;
    uint64_t hits = 0;
    for(auto it = prios.crbegin(); it != prios.crend(); it++) {
        if(curSize + it->size > cacheSize) {
            std::cout << cacheSize << " " << hits << " " << reqc << "\n";
            cacheSize*=2;
        }
        hits += it->reqCount - 1; // first request -> no prefetching
        curSize += it->size;
    }
    std::cout << -1 << " " << hits << " "  << reqc << "\n";

    return 0;
}



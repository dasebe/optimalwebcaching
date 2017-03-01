#include <cassert>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <functional>

typedef std::tuple<uint64_t,uint64_t,bool> traceEntry; //id, size, hasNext

int main(int argc, char* argv[]) {

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    
    // parse trace file
    std::vector<traceEntry> trace;
    std::map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0, uniqc=0;
    std::vector<std::vector<std::pair<uint64_t, uint64_t > > > sameTime;


    while(traceFile >> time >> id >> size) {
        if(lastSeen.count(std::make_pair(id,size))>0) {
            // see object second time
            std::get<2>(trace[lastSeen[std::make_pair(id,size)]]) = true;
            for(uint64_t i=lastSeen[std::make_pair(id,size)]; i<reqc; i++) {
                sameTime[i].push_back(std::make_pair(id,size));
            }
        } else {
            // see for first time
            uniqc++;
        }
        trace.emplace_back(id,size,false);
        sameTime.push_back( std::vector<std::pair<uint64_t, uint64_t>>() );
        lastSeen[std::make_pair(id,size)]=reqc++;
    }

    std::multimap<uint64_t, uint64_t > timeWidth;
    uint64_t thisTimeWidth;
    reqc = 0;
    for(auto & it: sameTime) {
        thisTimeWidth = 0;
        for(auto & vit: it) {
            std::cout << vit.first << "(" << vit.second << ") ";
            thisTimeWidth+=vit.second;
        }
        timeWidth.emplace(thisTimeWidth, reqc);
        std::cout << " TW " << thisTimeWidth << " " << reqc;
        std::cout << "\n";
        reqc++;
    }

    // largest width
    auto widthIt = --timeWidth.end();
    uint64_t deltaStar = widthIt->first;
    while (deltaStar > 0) {
        const uint64_t timeStar = widthIt->second;
        std::cout <<  deltaStar << " " <<  timeStar<< "\n";
        // sort request intervals by size
        std::sort(sameTime[timeStar].begin(), sameTime[timeStar].end(),
                  // lambda: sort using second pair entry (the size)
                  [](const std::pair<uint64_t,uint64_t> &left, const std::pair<uint64_t,uint64_t> &right) {
                return left.second < right.second;
            });
        for(auto & it : sameTime[timeStar]) {
            std::cout << it.first << "," << it.second << " | ";
        }
        std::cout << "\n";
        const uint64_t maxSize = (--sameTime[timeStar].end())->second;
        
        break;
        // again look at last/highest
        widthIt = --timeWidth.end();
        deltaStar = widthIt->first;
    }



    return 0;
}

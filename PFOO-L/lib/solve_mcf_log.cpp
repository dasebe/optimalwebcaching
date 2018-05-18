#include "solve_mcf.h"
#include <map>

void cacheAlg(std::vector<trEntry> & trace, uint64_t cacheSize) {
    std::unordered_map<std::pair<uint64_t, uint64_t>, trEntry*> cacheState;
    uint64_t currentSize = 0;

    for(size_t i=0; i<trace.size(); i++) {
        trEntry *cur = &trace[i];
        // check in cache
        if(cacheState.count(std::make_pair(cur->id,cur->size)) > 0) {
            // cache hit
            cur->hit = true;
        } else {
            // cache miss
            // admit if hasNext
            if(cur->hasNext && cur->size < cacheSize) {
                cacheState[std::make_pair(cur->id,cur->size)] = cur;
                currentSize += cur->size;
                LOG("admitted",cur->id,cur->size,currentSize);
                // evict if needed
                std::multimap<int64_t, trEntry *> distMap;
                int64_t curDistance;
                if(currentSize > cacheSize) {
                    for(auto & it: cacheState) {
                        trEntry * cand = it.second;
                        //LOG("cache Scan id",cand->id,cand->size,cand->nextSeen);
                        if(cand->nextSeen > i)
                            curDistance = cand->nextSeen - i;
                        else
                            curDistance = i - cand->nextSeen;
                        distMap.emplace(curDistance, cand);
                    }
                    while(currentSize > cacheSize) {
                        const auto victIt = --(distMap.end());
                        const trEntry * evictVictim = victIt->second;
                        LOG("vict",victIt->first,evictVictim->id,evictVictim->size);
                        cacheState.erase(std::make_pair(evictVictim->id,evictVictim->size));
                        currentSize -= evictVictim->size;
                        distMap.erase(victIt);
                    }
                }
            }
        }
    }
}

void printRes(std::vector<trEntry> & trace, std::string algName) {
    uint64_t hitc = 0, reqc = 0;
    for(auto & it: trace) {
        reqc++;
        if(it.hit)
            hitc++;
        LOG("tr",it.id,it.nextSeen,it.hit);
    }
    std::cout << algName << " hitc " << hitc << " reqc " << reqc << " ohr " << double(hitc)/reqc << std::endl;
}

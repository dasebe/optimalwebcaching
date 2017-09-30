#include <random>
#include "solve_mcf.h"

void cacheAlg(std::vector<trEntry> & trace, uint64_t cacheSize, size_t sampleSize) {
    std::unordered_map<std::pair<uint64_t, uint64_t>, trEntry*> cacheState;
    std::vector<size_t> cacheList;
    uint64_t currentSize = 0;

    std::default_random_engine generator;

    for(size_t i=0; i<trace.size(); i++) {
        trEntry *cur = &trace[i];
        // check in cache
        if(cacheState.count(std::make_pair(cur->id,cur->size)) > 0) {
            // cache hit
            cur->hit = true;
            cacheState[std::make_pair(cur->id,cur->size)] = cur;
        } else {
            // cache miss
            // admit if hasNext
            if(cur->hasNext && cur->size < cacheSize) {
                cacheState[std::make_pair(cur->id,cur->size)] = cur;
                cacheList.push_back(i);
                currentSize += cur->size;
                LOG("admitted",cur->id,cur->size,currentSize);
                
                // evict if needed
                while(currentSize > cacheSize) {
                    //                    LOG("eviction needed",currentSize,cacheSize,0);
                    // sampling-based eviction
                    std::uniform_int_distribution<size_t> distribution(0,cacheList.size()-1);

                    int64_t maxDistance = -1, curDistance; 
                    size_t victimIndex = 0;
                    for (size_t i=0; i<sampleSize; ++i) {
                        const size_t candIndex = distribution(generator);
                        trEntry * cand = &trace[cacheList[candIndex]];
                        //LOG("cache Scan id",cand->id,cand->size,cand->nextSeen);
                        if(cand->nextSeen > i)
                            curDistance = cand->nextSeen - i;
                        else
                            curDistance = i - cand->nextSeen;
                        if(curDistance > maxDistance) {
                            maxDistance = curDistance;
                            victimIndex = candIndex;
                        }
                    }
                    trEntry * evictVictim = &trace[cacheList[victimIndex]];
                    LOG("vict",maxDistance,evictVictim->id,evictVictim->size);
                    if(victimIndex != cacheList.size()-1) {
                        cacheList[victimIndex] = cacheList[cacheList.size()-1];
                    }
                    cacheList.pop_back();
                    if(cacheState.count(std::make_pair(evictVictim->id, evictVictim->size)) == 0) {
                        std::cerr << "BUG: in cacheList but not in cacheState " << victimIndex << "\n";
                    }
                    cacheState.erase(std::make_pair(evictVictim->id, evictVictim->size));
                    currentSize -= evictVictim->size;
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

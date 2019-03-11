#include <random>
#include "solve_mcf.h"

// belady's algorithm / min algorithm
void cacheAlg(std::vector<trEntry> & trace, uint64_t cacheSize, size_t sampleSize) {
    std::unordered_map<std::pair<uint64_t, uint64_t>, trEntry*> cacheState;
    std::vector<size_t> cacheList;
    uint64_t currentSize = 0;

    std::default_random_engine generator;

    for(size_t i=0; i<trace.size(); i++) {
        trEntry *cur = &trace[i];
        // check in cache
        if(cacheState.count(std::make_pair(cur->id,cur->size)) > 0) {
            // cache hit is fraction of bytes cached
            const auto oSize =cacheState[std::make_pair(cur->id,cur->size)]->size;
            const auto oCached = cacheState[std::make_pair(cur->id,cur->size)]->cachedBytes;
            if(oCached > oSize) {
                std::cerr << "PANIC " << oSize << " oC " << oCached << "\n";
            }
            cur->hit = oCached/double(oSize);
            cur->cachedBytes = oCached;
            // could increase cache fraction?
            if(currentSize < cacheSize && oSize > oCached) {
                // can increase fraction
                // std::cerr << "before hit : curCS " << currentSize << " CS " << cacheSize << " oS " << oSize << " oC " << oCached << "\n";
                if(cacheSize - currentSize > oSize - oCached) {
                    // cache full object
                    cur->cachedBytes = oSize;
                    currentSize += oSize - oCached;
                } else {
                    // increase fraction
                    cur->cachedBytes += cacheSize - currentSize;
                    currentSize = cacheSize;
                }
                // std::cerr << "after hit : curCS " << currentSize << " CS " << cacheSize << " oS " << oSize << " oC " << cur->cachedBytes << "\n";
            }

            // update to next object
            cacheState[std::make_pair(cur->id,cur->size)] = cur;
        } else {
            // cache miss
            // admit if hasNext
            if(cur->hasNext && cur->size < cacheSize && cur->size > 0) {
                // admit full object
                cacheState[std::make_pair(cur->id,cur->size)] = cur;
                cur->cachedBytes = cur->size;
                cacheList.push_back(i);
                currentSize += cur->size;
                LOG("admitted",cur->id,cur->size,currentSize);
                
                // evict if needed
                while(currentSize > cacheSize) {

                    // initialize with currently admitted
                    int64_t curDistance;
                    if(cur->nextSeen > i)
                        curDistance = (cur->nextSeen - i);
                    else
                        curDistance = (i - cur->nextSeen);
                    int64_t maxDistance = curDistance;
                    size_t victimIndex = cacheList.size()-1;
                    
                    // sample from rest
                    std::uniform_int_distribution<size_t> distribution(0,cacheList.size()-2);
                    for (size_t si=0; si<sampleSize; ++si) {
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
                    // evict a fraction if that's sufficient
                    if(currentSize - cacheSize >= evictVictim->cachedBytes) {
                        // evict full object
                        if(victimIndex != cacheList.size()-1) {
                            cacheList[victimIndex] = cacheList[cacheList.size()-1];
                        }
                        cacheList.pop_back();
                        if(cacheState.count(std::make_pair(evictVictim->id, evictVictim->size)) == 0) {
                            std::cerr << "BUG: in cacheList but not in cacheState " << victimIndex << "\n";
                        }
                        cacheState.erase(std::make_pair(evictVictim->id, evictVictim->size));
                        currentSize -= evictVictim->size;
                    } else {
                        //std::cerr << "before evict : curCS " << currentSize << " CS " << cacheSize << " oS " << evictVictim->size << " oC " << evictVictim->cachedBytes << "\n";
                        // evict partial object
                        evictVictim->cachedBytes -= currentSize - cacheSize;
                        currentSize = cacheSize;
                        //std::cerr << "after evict : curCS " << currentSize << " CS " << cacheSize << " oS " << evictVictim->size << " oC " << evictVictim->cachedBytes << "\n";
                    }
                }
            }
        }
    }
}

// just print out hit ratios
void printRes(std::vector<trEntry> & trace) {
    double hitc = 0, bytehitc = 0, byteSum = 0;
    for(auto & it: trace) {
        byteSum += it.size;

        // only count fraction
        hitc+=it.hit;
        bytehitc += it.hit * double(it.size);
        LOG("tr",it.id,it.nextSeen,it.hit);
    }
    std::cout << "BeladySplit ohr " << double(hitc)/trace.size() << " bhr " << double(bytehitc)/byteSum << "\n";
}

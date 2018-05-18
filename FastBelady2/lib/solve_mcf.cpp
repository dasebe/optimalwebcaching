#include "solve_mcf.h"
#include <cassert>
#include <cmath>

typedef std::multimap<long double, std::pair<uint64_t, uint64_t> > pMap_t;
typedef std::unordered_map<std::pair<uint64_t, uint64_t>, pMap_t::iterator> pRMap_t;

inline void checkInf(long double res) {
    if(std::isinf(res)) {
        std::cerr << "INF " << res << std::endl;
        assert(!std::isinf(res));
    }
}    


inline long double calcPrio(size_t i, trEntry* cand, long double baseL, long double normalizer) {
    long double res = 0;
    if(!cand->hasNext)
        res = normalizer;
    else if(cand->nextSeen > i)
        res = cand->nextSeen;
    else 
        res = -1.0L * i;
    checkInf(res);
    return res;
}


void cacheAlg(std::vector<trEntry*> & trace, uint64_t cacheSize) {
    std::unordered_map<std::pair<uint64_t, uint64_t>, trEntry*> cacheState;
    pMap_t prio;
    pRMap_t  revPrio;
    uint64_t currentSize = 0;

    long double baseL = 0;

    for(size_t i=0; i<trace.size(); i++) {
        LOG("-----start-----",i,trace.size(),0);
        trEntry *cur = trace[i];
        LOG("cur",i,trace.size(),0);
        const auto idSize = std::make_pair(cur->id,cur->size);
        // check in cache
        if(cacheState.count(idSize) > 0) {
            // cache hit
            LOG("hit",cur->id,cur->size,0);
            cur->hit = true;
            assert(revPrio.find(idSize) != revPrio.end());
            prio.erase(revPrio[idSize]);
            const long double curPrio = calcPrio(i, cur, baseL, trace.size());
            revPrio[idSize] = prio.emplace(curPrio, idSize);
            //            checkInf(revPrio[idSize]->first);
            LOG("hitted",cur->id,cur->size,revPrio[idSize]->first);
        } else {
            // cache miss
            // admit if hasNext
            if(cur->hasNext) {
                LOG("start admit",cur->id,cur->size,0);
                cacheState[idSize] = cur;
                currentSize += cur->size;
                const long double curPrio = calcPrio(i, cur, baseL, trace.size());
                revPrio[idSize] = prio.emplace(curPrio, idSize);
                //                checkInf(revPrio[idSize]->first);
                LOG("admitted",cur->id,cur->size,curPrio);
                // evict if needed
                while(currentSize > cacheSize) {
                    LOG("evict start",prio.size(),0,0);
                    pMap_t::iterator victimIt = --prio.end();
                    LOG("evict prio",victimIt->first,victimIt->second.first,victimIt->second.second);
                    //                    baseL += victimIt->first;
                    assert(cacheState.find(victimIt->second) != cacheState.end());
                    trEntry * evictVictim = cacheState[victimIt->second];
                    cacheState.erase(victimIt->second);
                    assert(revPrio.find(victimIt->second) != revPrio.end());
                    revPrio.erase(victimIt->second);
                    prio.erase(victimIt);
                    currentSize -= evictVictim->size;
                    LOG("evict end",currentSize,evictVictim->id,evictVictim->size);
                }
            }
        }
    }
}

void printRes(std::vector<trEntry*> & trace, std::string algName) {
    uint64_t hitc = 0, reqc = 0;
    for(auto it: trace) {
        reqc++;
        if(it->hit)
            hitc++;
        LOG("tr",it->id,it->nextSeen,it->hit);
    }
    std::cout << algName << " hitc " << hitc << " reqc " << reqc << " ohr " << double(hitc)/reqc << std::endl;
}

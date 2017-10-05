#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iomanip>
#include "lib/parse_trace.h"

using namespace lemon;


const bool feasibleCacheAll(const std::vector<trEntry > & trace, const std::vector<utilEntry> & utils, const int64_t cacheSize, const double minUtil) {

    // delta data structures: from localratio technique
    int64_t curDelta;

    curDelta = -cacheSize;
    std::unordered_set<size_t> indexCached;

    for(size_t i=0; i<utils.size(); i++) {
        const utilEntry & uEntry = utils[i];
        const trEntry & tEntry = trace[uEntry.index];

        // evict: if cached no next req
        if(( indexCached.count(uEntry.index)>0 ) && // cached
           (
            (uEntry.utility>=minUtil) || // not in util in interval
            (tEntry.nextSeen == 0)
            ) // no other request
           ) {
            LOG("evict",curDelta,uEntry.index,tEntry.nextSeen);
            curDelta -= tEntry.size;
        }

        // cache: if in utils bracket
        if(uEntry.utility > minUtil) {

            // if not cached
            if( indexCached.count(uEntry.index) == 0) {
                curDelta += tEntry.size;
                LOG("admit",curDelta,uEntry.index,tEntry.nextSeen);
            } // else: don't need update the size/width

            // add index to cache
            indexCached.insert(tEntry.nextSeen);

            // exit if we exceeded cache size
            if(curDelta > 0) {
                LOG("overflow",curDelta,uEntry.index,tEntry.nextSeen);
                return false;
            }
        }

        // delete current (now outdated) entry
        indexCached.erase(uEntry.index);
    }
    // all went well
    return true;

}

void initRemCap(const std::vector<trEntry > & trace, const std::vector<utilEntry> & utils, const int64_t cacheSize, const double minUtil, int64_t * remcap) {

    // delta data structures: from localratio technique
    int64_t curDelta;

    curDelta = -cacheSize;
    std::unordered_set<size_t> indexCached;

    for(size_t i=0; i<utils.size(); i++) {
        const utilEntry & uEntry = utils[i];
        const trEntry & tEntry = trace[uEntry.index];

        if(( indexCached.count(uEntry.index)>0 ) && // cached
           (
            (uEntry.utility>=minUtil) || // not in util in interval
            (tEntry.nextSeen == 0)
            ) // no other request
           ) {
            LOG("evict",curDelta,uEntry.index,tEntry.nextSeen);
            curDelta -= tEntry.size;
        }

        // cache: if in utils bracket
        if(uEntry.utility > minUtil) {

            // if not cached
            if( indexCached.count(uEntry.index) == 0) {
                curDelta += tEntry.size;
                LOG("admit",curDelta,uEntry.index,tEntry.nextSeen);
            } // else: don't need update the size/width

            // add index to cache
            indexCached.insert(tEntry.nextSeen);

            // exit if we exceeded cache size
            if(curDelta > 0) {
                std::cerr << "OVERFLOW\n";
            }
        }

        // delete current (now outdated) entry
        indexCached.erase(uEntry.index);

        // set remaining cap
        remcap[i] = -curDelta;
    }
}



int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << argv[0] << " traceFile cacheSize" << std::endl;
        return 1;
    }

    std::cerr << "starting..\n";

    std::string path(argv[1]);
    const int64_t cacheSize(std::stoull(argv[2]));

    // parse trace file
    std::vector<trEntry> trace;
    uint64_t uniqCount;
    uint64_t totalReqc = parseTraceFile(trace, path, cacheSize, uniqCount);
    std::cerr << "scanned trace n=" << totalReqc << std::endl;

    // ordered list of utilities and check that objects have size less than cache size

    std::vector<utilEntry> utils;
    std::cerr << "reserving " <<  trace.size() - uniqCount << " " << trace.size() << " " << uniqCount << "\n";
    utils.reserve(trace.size() - uniqCount);
    for(size_t i=0; i<trace.size(); i++) {
        const trEntry & cur = trace[i];
        if(cur.nextSeen > 0) {
            // calculate utility
            const double intervalLength = cur.nextSeen - i;
            const double utilityDenominator = cur.size*intervalLength;
            if(utilityDenominator<=0) {
                std::cerr << "fail " << cur.size << " ns " << cur.nextSeen << " " << i << " " << utilityDenominator << " " << intervalLength << "\n";
            }
            assert(utilityDenominator>0);
            utils.emplace_back(1.0L/utilityDenominator, i);
        }
    }
    std::sort(utils.begin(), utils.end());
    std::cerr << "ordered " << utils.size() << " utilities\n";

    // find threshold region: where some are not cacheable anymore
    size_t lOffset = 0, hOffset = 1e5, lowestUnfeasible = 0;
    bool foundInfeas = false;
    while(true) {
        std::cout << "searching thres " << lOffset << " " << utils[utils.size()-lOffset-1].utility << " " << hOffset  << " " << utils[utils.size()-hOffset].utility;
        const bool cacheAll = feasibleCacheAll(trace, utils, cacheSize, utils[utils.size()-hOffset].utility);
        std::cout << " - " << cacheAll  << "\n";
        if(cacheAll) {
            if(!foundInfeas) {
                lOffset = hOffset;
                hOffset *= 2;
                if(hOffset>utils.size()) {
                    hOffset = utils.size();
                }
            } else {
                lOffset = hOffset;
                const size_t stepSize = (lowestUnfeasible - lOffset)/2;
                if(stepSize < 1e3) {
                    break;
                }
                hOffset = lOffset + stepSize;
            }
        } else {
            lowestUnfeasible = hOffset;
            foundInfeas = true;
            const size_t stepSize = (hOffset - lOffset)/2;
            if(stepSize < 1e3) {
                break;
            }
            hOffset = lOffset + stepSize;
        }
    }

    //    std::cout << std::fixed << std::setprecision(20);
    std::cout << "found thres " << lOffset << " " << utils[utils.size()-lOffset].utility << " " << hOffset  << " " << utils[utils.size()-hOffset].utility << "\n";

    // create rem cap mask
    int64_t * remcap = new int64_t[trace.size()];
    initRemCap(trace, utils, cacheSize, utils[utils.size()-lOffset].utility, remcap);

    uint64_t hitc = lOffset, counter = 0;

    for (auto it = utils.rbegin()+(lOffset-1); it != utils.rend(); ++it) {
        
        const trEntry & cur = trace[it->index];
        bool enoughSpace = true;
        const auto s = cur.size;

        for(size_t i = it->index; i<cur.nextSeen; i++) {
            if(remcap[i] < s) {
                enoughSpace = false;
                break;
            }
        }

        if(enoughSpace) {
#pragma omp simd
            for(size_t i = it->index; i<cur.nextSeen; i++) {
                remcap[i] -= s;
            }
            hitc++;
        }

        if(counter++ > 100000) {
            std::cout << "lboundfast " << cacheSize << " " << hitc << " " <<  totalReqc << " " << (double)hitc/totalReqc << "\n";
            counter = 0;
        }
    }

            
    std::cout << "lboundfast " << cacheSize << " " << hitc << " " <<  totalReqc << " " << (double)hitc/totalReqc << "\n";

    return 0;
}

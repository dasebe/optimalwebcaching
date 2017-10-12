#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iomanip>
#include "lib/parse_trace.h"
#include <limits>

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << argv[0] << " traceFile" << std::endl;
        return 1;
    }

    std::cerr << "starting..\n";

    std::string path(argv[1]);

    // parse trace file
    std::vector<trEntry> trace;
    uint64_t uniqCount = 0;
    uint64_t totalReqc = parseTraceFile(trace, path, uniqCount);
    std::cerr << "scanned trace n=" << totalReqc << std::endl;

    // ordered list of utilities
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

    // create used-up capacity mask
    int64_t * usedcap = new int64_t[trace.size()];
#pragma omp simd
    for(size_t i=0; i<trace.size(); i++) {
        usedcap[i] = 0;
    }

    std::unordered_set<size_t> cached;
    uint64_t counter = 0;


    for(int64_t cs=1; cached.size() != trace.size()-uniqCount; cs = cs * 2) {

        // safe guard
        if(cs > std::numeric_limits<int64_t>::max()/2) {
            break;
        }

        std::cerr << cs;
        for (auto it = utils.rbegin(); it != utils.rend(); it++) {

            if(cached.count(it->index)>0) {
                continue; // already cached
            }
            
            const trEntry & cur = trace[it->index];
            const auto s = cur.size;
            if(s > cs) {
                continue; // skip if too large to store at the moment
            }

            bool cacheThis = true;
            for(size_t i = it->index; i<cur.nextSeen; i++) {
                if(usedcap[i] + s > cs) {
                    cacheThis = false;
                    break;
                }
            }

            if(cacheThis) {
#pragma omp simd
                for(size_t i = it->index; i<cur.nextSeen; i++) {
                    usedcap[i] += s;
                }
                cached.insert(it->index);
            }

            if(counter++ > 10000) {
                std::cerr << ".";
                counter = 0;
            }

        }


        std::cerr << std::endl;
        counter = 0;
        std::cout << "lboundfastc " + std::to_string(cs) + " " + std::to_string(cached.size()) + " " + std::to_string(totalReqc) +  " " + std::to_string((double)(cached.size())/totalReqc) + "\n";
    }

    return 0;
}

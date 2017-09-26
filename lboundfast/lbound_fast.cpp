#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "lib/parse_trace.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << argv[0] << " traceFile cacheSize" << std::endl;
        return 1;
    }

    std::cerr << "starting..\n";

    std::string path(argv[1]);
    int64_t cacheSize(std::stoull(argv[2]));

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
            assert(utilityDenominator>0);
            utils.emplace_back(1.0L/utilityDenominator, i);
        }
    }
    std::sort(utils.begin(), utils.end());
    std::cerr << "ordered " << utils.size() << " utilities\n";

    int64_t * remcap = new int64_t[trace.size()];
    std::fill(remcap, remcap + trace.size(), cacheSize);
    uint64_t hitc = 0, counter = 0;

    for (auto it = utils.rbegin(); it != utils.rend(); ++it) {
        //        std::cerr << it->utility << " " << it->index << "\n";
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
            for(size_t i = it->index; i<cur.nextSeen; i++) {
                remcap[i] -= s;
            }
            hitc++;
        }

        if(counter++ > 1000000) {
            std::cerr << "step " << hitc << "\n";
            counter = 0;
        }
    }

            
    std::cout << "lboundfast " << cacheSize << " " << hitc << " " <<  totalReqc << " " << (double)hitc/totalReqc << "\n";

    return 0;
}

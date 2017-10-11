#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <set>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 5) {
        std::cerr << argv[0] << " traceFile cacheSize solverParam maxEjectSize" << std::endl;
        return 1;
    }

    std::cerr << "starting..\n";

    std::string path(argv[1]);
    uint64_t cacheSize(std::stoull(argv[2]));
    int solverPar(std::stoi(argv[3]));
    uint64_t maxEjectSize(std::stoull(argv[4]));

    // parse trace file
    std::vector<trEntry> trace;
    uint64_t totalUniqC = parseTraceFile(trace, path, cacheSize);
    uint64_t totalReqc = trace.size();
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc - totalUniqC) {
        maxEjectSize = totalReqc - totalUniqC;
    }
    
    // ordered list of utilities and check that objects have size less than cache size
    std::multiset<double> utilities;
    for(auto & it: trace) {
        if(it.nextSeen > 0) {
            assert(it.utility>=0);
            utilities.insert(it.utility);
        }
    }

    // get utility boundaries for ejection sets (based on ejection set size)
    std::vector<double> utilSteps;
    utilSteps.push_back(1); // max util as start
    uint64_t curEjectSize = 0;
    LOG("ejSize",maxEjectSize,trace.size(),utilities.size());
    assert(maxEjectSize>0);
    for(auto it=--(utilities.end()); it!=utilities.begin(); it--) { //TBD last entry?
        curEjectSize++;
        if(curEjectSize >= maxEjectSize/2 && (*it != *(--(utilSteps.end())) ) ) {
            utilSteps.push_back(*it);
            //DEBUG
            LOG("utilStep",*it,0,curEjectSize);
            curEjectSize = 0;
        }
    }
    utilSteps.push_back(0); // min util as end
    // clear out utilities
    std::multiset<double>().swap(utilities);
    std::cerr << "ejection sets - #sets: " << utilSteps.size() << " |set|: " << maxEjectSize << "\n";
        

    long double curCost=0, curHits, overallHits = 0;
    uint64_t integerHits = 0;
    size_t effectiveEjectSize=0;
    
    // binary search for highest k such that (k+1) unfeasible and (k) feasible to set all dvars with utility > minUtil[k] to 1
    size_t kMinUnFeasible = utilSteps.size()-3;
    size_t kMaxFeasible = 0;
    bool soFarFeasibleCacheAll = true;

    std::cout.precision(10);
    std::cout << std::fixed;

    // LNS iteration steps
    for(size_t k=kMaxFeasible; k+2<utilSteps.size(); k++) {
        // set step's util boundaries
        const double minUtil = utilSteps[k+2];
        const double maxUtil = utilSteps[k];

        std::unordered_map<size_t, int> traceToArc;

        // check if we can simply add all intervals of the current ejection set
        if(soFarFeasibleCacheAll) {
            soFarFeasibleCacheAll = feasibleCacheAll(trace, cacheSize, minUtil);

            // output iteration statistics
            std::cerr << "k " << k << " lU " << minUtil << " uU " << 1
                      << " feas " << soFarFeasibleCacheAll << " maxfeas " << kMaxFeasible << " minunfeas " << kMinUnFeasible
                      << " pp " << 1 << " oR " << totalReqc << "\n";
                
            // not feasible
            if(!soFarFeasibleCacheAll) {
                // update kMinUnFeasible as this was not feasible
                if(k<kMinUnFeasible) {
                    kMinUnFeasible = k;
                }
                // update k to middle of interval
                k = kMaxFeasible + (kMinUnFeasible-kMaxFeasible)/2-1;
                // give it another try
                soFarFeasibleCacheAll = true;
            } else { // feasible
                // update kMaxFeasible as this was feasible
                if(k>kMaxFeasible)
                    kMaxFeasible = k;
                // update k to middle of interval
                k = kMaxFeasible + (kMinUnFeasible-kMaxFeasible)/2-1;
            }
            // last iteration
            if(kMaxFeasible + 1 >= kMinUnFeasible) {
                // set all dvars of injection (so far) to 1
                curCost = 0;
                curHits = 0;
                overallHits = 0;
		integerHits = 0;
                effectiveEjectSize = 0;
                // set k to max known feasible sol
                k=kMaxFeasible;
                const double minUtil = utilSteps[k];
                for(uint64_t i=0; i<trace.size(); i++) {
                    effectiveEjectSize++;
                    if(trace[i].utility>=minUtil && cacheSize >= trace[i].size) {
                        trace[i].dvar = 1;
                        curHits++;
                        overallHits++;
			integerHits++;
                    }
                }
                // output iteration statistics
                std::cout << "k " << k << " lU " << minUtil << " uU " << 1
                          << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                          << " oH " << overallHits << " oR " << totalReqc << " iH " << integerHits << "\n";
                // we now continue with MCF
                soFarFeasibleCacheAll = false;
            }
            // we can now skip the rest of this iteration
            continue;
        }

        // create MCF digraph with arc utilities in [minUtil,maxUtil]
        SmartDigraph g; // mcf graph
        SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
        SmartDigraph::ArcMap<double> cost(g); // mcf costs
        SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies
        effectiveEjectSize = createMCF(g, trace, cacheSize, cap, cost, supplies, minUtil, maxUtil, traceToArc);

        // solve this MCF
        SmartDigraph::ArcMap<uint64_t> flow(g);
        curCost = solveMCF(g, cap, cost, supplies, flow, solverPar);

        // write DVAR to trace
        curHits = 0;
        //        overallHits = 0;
        const double middleUtil = utilSteps[k+1];
        for(auto it: traceToArc) {
            const auto arcId = it.second;
            if(trace[it.first].utility>=middleUtil && trace[it.first].utility<maxUtil && cacheSize >= trace[it.first].size) {
                trace[it.first].dvar = 1.0L - flow[g.arcFromId(arcId)]/static_cast<long double>(trace[it.first].size);
                curHits += trace[it.first].dvar;
                LOG("dv",it.first,trace[it.first].dvar,trace[it.first].size);
                assert(trace[it.first].dvar >= 0 && trace[it.first].dvar<=1);
                overallHits += trace[it.first].dvar;
		if(trace[it.first].dvar > 0.99) {
		  integerHits++;
		}
            }
        }

        // output iteration statistics
        std::cout << "k " << k << " lU " << minUtil << " uU " << maxUtil
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
		  << " oH " << overallHits << " oR " << totalReqc << " iH " << integerHits << "\n";
    }

    /// final summing up
    overallHits = 0;
    for(uint64_t i=0; i<trace.size(); i++) {
        overallHits += trace[i].dvar;
    }
    // output iteration statistics
    std::cout << "final " << -1 << " lU " << -1 << " uU " << -1
              << " cC " << curCost << " cH " << -1 << " cR " << -1
              << " oH " << overallHits << " oR " << totalReqc << "\n";

    return 0;
}




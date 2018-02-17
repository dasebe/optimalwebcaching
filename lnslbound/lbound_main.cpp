#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <set>
#include <iomanip>
#include <algorithm>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 6) {
        std::cerr << argv[0] << " traceFile cacheSize solverParam maxEjectSize resultPath" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    uint64_t cacheSize(std::stoull(argv[2]));
    int solverPar(std::stoi(argv[3]));
    uint64_t maxEjectSize(std::stoull(argv[4]));
    std::string resultPath(argv[5]);

    // parse trace file
    std::vector<trEntry> trace;
    uint64_t totalUniqC = parseTraceFile(trace, path);
    uint64_t totalReqc = trace.size();
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << " vs " << trace.size() << std::endl;

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc - totalUniqC) {
        maxEjectSize = totalReqc - totalUniqC;
    }
    
    // ordered list of utilities and check that objects have size less than cache size
    std::vector<double> utilSteps2;
    for(auto & it: trace) {
        if(it.size > cacheSize) {
            it.hasNext = false;
        }
        if(it.hasNext) {
            assert(it.utility>=0);
            utilSteps2.push_back(it.utility);
        }
    }
    std::cerr << "scanned utilities " << utilSteps2.size() << std::endl;

    std::sort(utilSteps2.begin(), utilSteps2.end(),std::greater<double>());
    std::cerr << "sorted utilities " << utilSteps2.size() << std::endl;

    // get utility boundaries for ejection sets (based on ejection set size)
    std::vector<double> utilSteps;
    utilSteps.push_back(1); // max util as start
    uint64_t curEjectSize = 0;
    LOG("ejSize",maxEjectSize,trace.size(),utilities.size());
    assert(maxEjectSize>0);
    for(auto & it: utilSteps2) {
        curEjectSize++;
        if(curEjectSize >= maxEjectSize/2 && (it != *(--(utilSteps.end())) ) ) {
            utilSteps.push_back(it);
            //DEBUG
            LOG("utilStep",it,0,curEjectSize);
            curEjectSize = 0;
        }
    }
    utilSteps.push_back(0); // min util as end
    utilSteps2.clear();
    utilSteps2.shrink_to_fit();
    std::cerr << "defined ejection sets - #sets: " << utilSteps.size() << " |set|: " << maxEjectSize << "\n";

    long double curCost=0, curHits, overallHits;
    uint64_t integerHits = 0;
    size_t effectiveEjectSize=0;
    
    // LNS iteration steps
    for(size_t k=0; k+2<utilSteps.size(); k++) {

        // set step's util boundaries
        const double minUtil = utilSteps[k+2];
        const double maxUtil = utilSteps[k];

        std::cerr << "k1. " << k << " lU " << minUtil << " uU " << maxUtil
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                  << " oH " << overallHits << " oR " << totalReqc  << " iH " << integerHits << std::endl;


        // create MCF digraph with arc utilities in [minUtil,maxUtil]
        SmartDigraph g; // mcf graph
        SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
        SmartDigraph::ArcMap<double> cost(g); // mcf costs
        SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies
        effectiveEjectSize = createMCF(g, trace, cacheSize, cap, cost, supplies, minUtil, maxUtil);

        std::cerr << "k2. " << k << " lU " << minUtil << " uU " << maxUtil
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                  << " oH " << overallHits << " oR " << totalReqc  << " iH " << integerHits << std::endl;


        // solve this MCF
        SmartDigraph::ArcMap<uint64_t> flow(g);
        curCost = solveMCF(g, cap, cost, supplies, flow, solverPar);

        std::cerr << "k3. " << k << " lU " << minUtil << " uU " << maxUtil
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                  << " oH " << overallHits << " oR " << totalReqc  << " iH " << integerHits << std::endl;


        // write DVAR to trace
        curHits = 0;
        overallHits = 0;
        integerHits = 0;
        for(uint64_t i=0; i<trace.size(); i++) {
            if(trace[i].active) {
                trace[i].dvar = 1.0L - flow[g.arcFromId(trace[i].arcId)]/static_cast<long double>(trace[i].size);
                trace[trace[i].nextSeen].hit = trace[i].dvar;
                curHits += trace[i].dvar;
            }
            LOG("dv",i,trace[i].dvar,trace[i].size);
            assert(trace[i].dvar >= 0 && trace[i].dvar<=1);
            overallHits += trace[i].dvar;
            if(trace[i].dvar > 0.99) {
                integerHits++;
            }
        }

        // output iteration statistics
        std::cout << "k " << k << " lU " << minUtil << " uU " << maxUtil
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                  << " oH " << std::setprecision(20) << overallHits << " oR " << totalReqc  << " iH " << integerHits << std::endl;
    }

    // output decision variables and utilities
    std::ofstream resultfile(resultPath);

    for(auto & it: trace) {
        resultfile 
                   << it.id << " " << it.size << " "
                   << it.utility << " "
                   << it.dvar << " "
                   << it.hit << std::endl;
    }


    
    return 0;
}



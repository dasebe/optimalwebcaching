#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <set>
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
    uint64_t curHits = 0;
    uint64_t integerHits = 0;
    uint64_t floatHits = 0;
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // LNS iteration steps

    for(size_t k=0; k<trace.size(); k+=maxEjectSize/2) {

        // create MCF digraph with arc utilities in [minUtil,maxUtil]
        SmartDigraph g; // mcf graph
        SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
        SmartDigraph::ArcMap<double> cost(g); // mcf costs
        SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies
        const size_t effectiveEjectSize = createMCF(g, trace, cacheSize, cap, cost, supplies, k, k+maxEjectSize);

        // solve this MCF
        SmartDigraph::ArcMap<uint64_t> flow(g);
        const long double curCost = solveMCF(g, cap, cost, supplies, flow, solverPar);

        // write DVAR to trace
        curHits = 0;
        integerHits = 0;
        floatHits = 0;
        for(uint64_t i=0; i<trace.size(); i++) {
            if(trace[i].active) {
                trace[i].dvar = 1.0L - flow[g.arcFromId(trace[i].arcId)]/static_cast<long double>(trace[i].size);
                trace[trace[i].nextSeen].hit = trace[i].dvar;
                curHits += trace[i].dvar;
            }
            LOG("dv",i,trace[i].dvar,trace[i].size);
            assert(trace[i].dvar >= 0 && trace[i].dvar<=1);
            floatHits += trace[i].dvar;
            if(trace[i].dvar > 0.99)
                integerHits ++;
        }

        // output iteration statistics
        std::cout << "k " << k << " km " << k+maxEjectSize/2 << " uU " << 0
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                  << " fH " << floatHits << " oR " << totalReqc << " iH " << integerHits << "\n";
    }

    // output decision variables and utilities
    std::ofstream resultfile(resultPath);

    for(auto & it: trace) {
        resultfile 
                   << it.id << " " << it.size << " "
                   << it.dvar << " "
                   << it.hit << "\n";
    }


    // write DVAR to trace
    curHits = 0;
    integerHits = 0;
    floatHits = 0;
    for(uint64_t i=0; i<trace.size(); i++) {
        floatHits += trace[i].dvar;
        if(trace[i].dvar > 0.99)
            integerHits ++;
    }

    // output iteration statistics
    std::cout << "final " << -1 << " check " << trace.size() << " uU " << 0
              << " cC " << 0 << " cH " << curHits << " cR " << 0
              << " fH " << floatHits << " oR " << totalReqc << " iH " << integerHits << "\n";

    
    return 0;
}



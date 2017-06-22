#include <fstream>
#include <lemon/lgf_writer.h>
#include <cassert>
#include <unordered_set>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << argv[0] << " traceFile cacheSize maxEjectSize" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    uint64_t maxEjectSize(std::stoull(argv[3]));
    //    std::string resultPath(argv[4]);

    // parse trace file
    std::vector<trEntry> trace;
    uint64_t totalUniqC = parseTraceFile(trace, path);
    uint64_t totalReqc = trace.size();
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // create mcf instance
    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    SmartDigraph::Node lastNode = createMCF(g, trace, cacheSize, cap, cost, supplies);
    
    std::cerr << "created graph with ";

    // dual-feasible solution: all alpha and pi are 0
    uint64_t nodes=0, vertices=0;
    double dualSol = 0;
    SmartDigraph::NodeMap<double> pi(g);
    SmartDigraph::ArcMap<double> alpha(g);

    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodes;
        pi[n] = 0;
        dualSol += pi[n] * supplies[n];
    }
    std::cerr << nodes << " nodes ";
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++vertices;
        alpha[a] = 0;
        dualSol += alpha[a] * cap[a];
    }
    std::cerr << vertices << " arcs with dual cost " << dualSol << std::endl;

    for (SmartDigraph::OutArcIt a(g, lastNode); a!=INVALID; ++a) ++nodes;

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc - totalUniqC) {
        maxEjectSize = totalReqc - totalUniqC;
    }
    
    // LNS iteration steps
    std::unordered_set<int> outerArcs;
    std::unordered_set<int> innerArcs;
    std::unordered_set<int> lArcs;
    std::unordered_set<int> rArcs;
    for(size_t kmin=0; kmin<trace.size(); kmin+=maxEjectSize/2) {
        const size_t kmax = kmin+maxEjectSize;
        outerArcs.clear();
        innerArcs.clear();
        lArcs.clear();
        rArcs.clear();
        for(uint64_t i=kmin; i<kmax; i++) {
            trEntry & curEntry = trace[i];
            outerArcs.insert(curEntry.outerArcId);
            outerArcs.insert(curEntry.innerArcId);
            SmartDigraph::Node lNode = g.source(g.arcFromId(curEntry.innerArcId));
            for (SmartDigraph::OutArcIt a(g, lNode); a!=INVALID; ++a) {
                const int curArcId = g.id(a);
                if(outerArcs.count(curArcId) == 0 && innerArcs.count(curArcId) == 0) {
                    lArcs.insert(curArcId);
                }
            }
            SmartDigraph::Node rNode = g.target(g.arcFromId(curEntry.innerArcId));
            for (SmartDigraph::OutArcIt a(g, rNode); a!=INVALID; ++a) {
                const int curArcId = g.id(a);
                if(outerArcs.count(curArcId) == 0 && innerArcs.count(curArcId) == 0) {
                    rArcs.insert(curArcId);
                }
            }
        }
    }


    return 0;
}

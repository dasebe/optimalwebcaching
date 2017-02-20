#include <lemon/lgf_writer.h>
#include <cassert>
#include "parse_trace.h"
#include "solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << argv[0] << " traceFile cacheSize solverParam" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    int solverPar(atoi(argv[3]));

    // parse trace file
    std::vector<traceEntry> trace;
    uint64_t totalUniqC = parseTraceFile(trace, path);
    uint64_t totalReqc = trace.size();
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // create mcf instance
    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    createMCF(g, trace, cacheSize, cap, cost, supplies);
    
    std::cerr << "created graph with ";
    uint64_t nodes=0, vertices=0;
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) ++nodes;
    std::cerr << nodes << " nodes ";
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) ++vertices;
    std::cerr << vertices << " arcs " << std::endl;

    SmartDigraph::ArcMap<uint64_t> flow(g);
    double solval = solveMCF(g, cap, cost, supplies, flow, solverPar);
    assert(solval>0);

    std::cerr << "solution par " << solverPar << " cost " << solval << " teqs " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << std::endl;
    
    for(auto & it: trace) {
        const uint64_t id=std::get<0>(it);
        const uint64_t size=std::get<1>(it);
        const uint64_t time=std::get<3>(it);
        const int arcId=std::get<4>(it);
        std::cout << time << " " << id << " " << size << " ";
        if(arcId==-1) 
            std::cout << "0\n";
        else
            std::cout << (size-flow[g.arcFromId(arcId)])/static_cast<double>(size) << "\n";
    }


    return 0;
}

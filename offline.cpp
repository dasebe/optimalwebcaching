#include <lemon/lgf_writer.h>
#include <lemon/network_simplex.h>
#include <lemon/capacity_scaling.h>
#include "offline.h"

using namespace lemon;

#ifdef NETWORKSIMPLEX
typedef NetworkSimplex<SmartDigraph, int64_t, double> SolverType;
#else
typedef CapacityScaling<SmartDigraph, int64_t, double, CapacityScalingDefaultTraits<SmartDigraph, int64_t, double> > SolverType;
#endif

// missing scale etc

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << argv[0] << " traceFile cacheSize CASscale" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    int scale(atoi(argv[3]));

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

    // solve the mcf instance
    SolverType solver(g);
    solver.upperMap(cap).costMap(cost).supplyMap(supplies);

    SolverType::ProblemType res;
#ifdef NETWORKSIMPLEX
    switch(scale) {
        case 1: res = solver.run(SolverType::FIRST_ELIGIBLE);
            std::cerr << "solver: NetworkSimplex FE ";
            break;
        case 2: res = solver.run(SolverType::BEST_ELIGIBLE);
            std::cerr << "solver: NetworkSimplex BE ";
            break;
        case 4: res = solver.run(SolverType::CANDIDATE_LIST);
            std::cerr << "solver: NetworkSimplex CL ";
            break;
        default: res = solver.run(SolverType::BLOCK_SEARCH);
            std::cerr << "solver: NetworkSimplex BS ";
            break;
    }
#else
    res = solver.run(scale);
    std::cerr << "solver: CapacityScaling S" << scale << " ";
#endif


    if(res==SolverType::INFEASIBLE) {
        std::cerr << "infeasible mcf" << std::endl;
        return -1;
    } else if (res==SolverType::UNBOUNDED) {
        std::cerr << "unbounded mcf" << std::endl;
        return -1;
    }

    std::cerr << "optimal solution: cost " << solver.totalCost<double>() << " teqs " << totalReqc << " OHR " << 1.0-(static_cast<double>(solver.totalCost<double>())+totalUniqC)/totalReqc << std::endl;
    
    SmartDigraph::ArcMap<uint64_t> flow(g);
    solver.flowMap(flow);
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

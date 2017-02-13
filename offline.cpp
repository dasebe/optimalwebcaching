#include <lemon/lgf_writer.h>
#include <lemon/capacity_scaling.h>
#include "offline.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));

    // parse trace file
    std::vector<std::tuple<uint64_t,uint64_t,bool> > trace;
    parseTraceFile(trace, path);

    // create mcf instance
    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int> supplies(g); // mcf demands/supplies
    createMCF(g, trace, cacheSize, cap, cost, supplies);
    trace.clear();
    
    std::cout << "created graph with ";

    int nodes=0, vertices=0;
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) ++nodes;
    std::cout << nodes << " nodes ";
    
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) ++vertices;
    std::cout << vertices << " arcs ";

    digraphWriter(g, std::cout).
        arcMap("capacity", cap).       // write cap into 'capacity'
        arcMap("cost", cost).          // write 'cost' for for arcs
        run();
    
    return 0;
}

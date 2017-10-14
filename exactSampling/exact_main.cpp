#include <fstream>
#include <lemon/lgf_writer.h>
#include <cassert>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 6) {
        std::cerr << argv[0] << " traceFile cacheSize solverParam resultPath samples" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    int solverPar(atoi(argv[3]));
    std::string resultPath(argv[4]);
    uint64_t samples(atoll(argv[5]));

    // parse trace file
    std::vector<traceEntry> trace;
    uint64_t totalUniqC = parseTraceFile(trace, path, samples);
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
    //    assert(solval>0);

   
    std::ofstream resultfile(resultPath);

    long double floatHits = 0;
    uint64_t integerHits = 0;

    for(auto & it: trace) {
        const uint64_t id=std::get<0>(it);
        const uint64_t size=std::get<1>(it);
        const uint64_t time=std::get<3>(it);
        const int arcId=std::get<4>(it);
        resultfile << time << " " << id << " " << size << " ";
        if(arcId==-1) 
            resultfile << "0\n";
        else {
            const long double dvar = (size-flow[g.arcFromId(arcId)])/static_cast<double>(size);
            resultfile << dvar << "\n";
            floatHits += dvar;
            if(dvar > 0.99) {
                integerHits++;
            }
        }
                                  
    }

    std::cout.precision(12);
    std::cout << std::fixed;

    std::cerr << "ExLP" << solverPar << " " << cacheSize << " hitc " << totalReqc-totalUniqC-solval << " reqc " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << " " << floatHits << " " << integerHits << std::endl;
    std::cout << "ExLP" << solverPar << " " << cacheSize << " hitc " << totalReqc-totalUniqC-solval << " reqc " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << " " << floatHits << " " << integerHits << std::endl;


    return 0;
}

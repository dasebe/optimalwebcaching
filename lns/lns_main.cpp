#include <lemon/lgf_writer.h>
#include <cassert>
#include <set>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

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
    std::vector<trEntry> trace;
    uint64_t totalUniqC = parseTraceFile(trace, path);
    uint64_t totalReqc = trace.size();
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // create mcf instance
    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    std::set<long double> utilities;
    long double minU=1, maxU=0;
    for(auto & it: trace) {
        utilities.insert(it.utility);
        minU=it.utility < minU ? it.utility : minU;
        maxU=it.utility > maxU ? it.utility : maxU;
    }
    std::cout << "minU " << minU << " maxU " << maxU << "\n";
    auto curMinUtility = utilities.begin();
    uint64_t maxEjectSize = trace.size()/10;
    uint64_t curEjectSize;
    bool lastIteration = false;

    for(uint64_t i=0; i<trace.size(); i++ && !lastIteration) {
        curEjectSize = 0;
        while(curEjectSize < maxEjectSize) {
            curMinUtility++;
            curEjectSize++;
            if(curMinUtility==utilities.end()) {
                lastIteration = true;
                break;
            }
        }
        std::cout << i << " " << *curMinUtility <<"\n";
    }
    
    cacheSize++;
    solverPar++;

    return 0;
}




        std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, int> > lastSeen;
        SmartDigraph::Arc curArc;
        SmartDigraph::Node curNode = g.addNode(); // initial node
        SmartDigraph::Node prevNode;

        // iterate over trace
        for(uint64_t i=0; i<trace.size(); i++) {
            const trEntry & curEntry = trace[i];
            // first: check if previous interval ended here
            if(lastSeen.count(std::make_pair(curEntry.id,curEntry.size))>0) {
                // create "outer" request arc
                const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[std::make_pair(curEntry.id,curEntry.size)].second);
                curArc = g.addArc(lastReq,curNode);
                cap[curArc] = curEntry.size;
                cost[curArc] = 1/static_cast <double>(curEntry.size);
                supplies[lastReq] += curEntry.size;
                supplies[curNode] -= curEntry.size;
                trace[lastSeen[std::make_pair(curEntry.id,curEntry.size)].first].arcId = g.id(curArc);
            }
            // second: if there is another request for this object
            if(curEntry.hasNext) {
                // save prev node as anchor for future arcs
                prevNode = curNode;
                lastSeen[std::make_pair(curEntry.id,curEntry.size)]=std::make_pair(i,g.id(prevNode));
                // create another node, "inner" capacity arc
                curNode = g.addNode(); // next node
                curArc = g.addArc(prevNode,curNode);
                cap[curArc] = cacheSize; 
                cost[curArc] = 0;
            }
        }

        SmartDigraph::ArcMap<uint64_t> flow(g);
        double solval = solveMCF(g, cap, cost, supplies, flow, solverPar);
        //    assert(solval>0);

        std::cerr << "solution par" << " " << solverPar << " cost " << solval << " teqs " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << std::endl;

        for(auto & it: trace) {
            std::cout << it.origTime << " " << it.id << " " << it.size << " ";
            if(it.arcId==-1) 
                std::cout << "0\n";
            else
                std::cout << (it.size-flow[g.arcFromId(it.arcId)])/static_cast<double>(it.size) << "\n";
        }


    }


#include <lemon/lgf_writer.h>
#include <cassert>
#include <vector>
#include <set>
#include <cmath>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"


// uncomment to enable debugging:
//#define DEBUG 1

#ifdef DEBUG
#define LOG(m,x,y,z) log_message(m,x,y,z,"\n")
#else
#define LOG(m,x,y,z)
#endif
inline void log_message(std::string m, double x, double y, double z, std::string e) {
    std::cerr << m << "," << x << "," << y  << "," << z << e;
}


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
    std::cerr << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc - totalUniqC) {
        maxEjectSize = totalReqc - totalUniqC;
    }
    
    long double curCost, curHits, overallHits;
    size_t effectiveEjectSize;
    std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, int> > lastSeen;
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curNode;
    SmartDigraph::Node prevNode;
    
    // LNS iteration steps
    size_t k=0;
    for(size_t kmin=0; kmin<trace.size(); kmin+=maxEjectSize/2) {
        // set step's boundaries
        const long double kmax = kmin+maxEjectSize;

        // create a graph with arcs between kmin,kmax
        // mcf instance data
        SmartDigraph g; // mcf graph
        SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
        SmartDigraph::ArcMap<double> cost(g); // mcf costs
        SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

        curNode = g.addNode(); // first node in interval

        effectiveEjectSize = 0;

        for(uint64_t i=0; i<trace.size(); i++) {
            trEntry & curEntry = trace[i];
            curEntry.active = false;

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
                trace[lastSeen[std::make_pair(curEntry.id,curEntry.size)].first].active = true;
                effectiveEjectSize++;
                // DEBUG
                LOG("oA",i,lastSeen[std::make_pair(curEntry.id,curEntry.size)].first,g.id(curArc));
                // delete lastSeen entry, as the next interval might not be part
                lastSeen.erase(std::make_pair(curEntry.id,curEntry.size));
            } 

            // in injection set IF hasnext
            if(curEntry.hasNext) {
                // AND EITHER interval starts before and ends after kmin
                if (i<kmin && curEntry.nextSeen > kmin) {
                    // save curNode (= first node) as anchnor for future arcs
                    lastSeen[std::make_pair(curEntry.id,curEntry.size)]=std::make_pair(i,g.id(curNode));
                } else
                    // OR interval starts in [kmin,kmax]
                    if( i>=kmin && i<kmax ) {
                        // save prev node as anchor for future arcs
                        prevNode = curNode;
                        lastSeen[std::make_pair(curEntry.id,curEntry.size)]=std::make_pair(i,g.id(prevNode));
                        // create another node, "inner" capacity arc
                        curNode = g.addNode(); // next node
                        curArc = g.addArc(prevNode,curNode);
                        cap[curArc] = cacheSize; 
                        // DEBUG
                        LOG("iA",i,cap[curArc],0);
                        cost[curArc] = 0;
                    }
            }
        }

        lastSeen.clear();
        std::cerr << "k " << k << " kmin " << kmin << " kmax "
                  << " starting 0 solver 0 x 0 cR " << effectiveEjectSize << "\n";

        
        // solve this MCF
        SmartDigraph::ArcMap<uint64_t> flow(g);
        curCost = solveMCF(g, cap, cost, supplies, flow, solverPar);

        // write DVAR to trace
        curHits = 0;
        overallHits = 0;
        for(uint64_t i=0; i<trace.size(); i++) {
            if(trace[i].active) {
                trace[i].dvar = 1.0L - flow[g.arcFromId(trace[i].arcId)]/static_cast<long double>(trace[i].size);
                curHits += trace[i].dvar;
            }
            LOG("dv",i,trace[i].dvar,trace[i].size);
            assert(trace[i].dvar >= 0 && trace[i].dvar<=1);
            overallHits += trace[i].dvar;
        }
        
        // output iteration statistics
        std::cout << "k " << k++ << " kmin " << kmin << " kmax " << kmax
                  << " cC " << curCost << " cH " << curHits << " cR " << effectiveEjectSize
                  << " oH " << overallHits << " oR " << totalReqc << "\n";
    }

    // output decision variables and utilities
    std::ofstream resultfile(resultPath);

    for(auto & it: trace) {
        resultfile << it.origTime << " "
                   << it.id << " " << it.size << " "
                   << it.nextSeen << " "
                   << it.dvar << "\n";
    }


    
    return 0;
}



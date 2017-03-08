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
    std::cout << m << "," << x << "," << y  << "," << z << e;
}

inline bool isInEjectSet(const long double minUtil, const long double maxUtil, const trEntry & curEntry) {
    return( curEntry.utility>=minUtil && curEntry.utility<maxUtil );
}


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
    std::cout << "scanned trace n=" << totalReqc << " m=" << totalUniqC << std::endl;

    // create mcf instance
    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    // ordered list of utilities and check that objects have size less than cache size
    std::multiset<long double> utilities;
    for(auto & it: trace) {
        if(it.size > cacheSize) {
            it.hasNext = false;
        }
        if(it.hasNext) {
            utilities.insert(it.utility);
        }
    }

    // get utility boundaries for ejection sets (based on ejection set size)
    std::vector<long double> utilSteps;
    utilSteps.push_back(1); // max util as start
    uint64_t maxEjectSize = utilities.size()/std::sqrt(totalReqc-totalUniqC); //+1
    uint64_t curEjectSize = 0;
    LOG("ejSize",maxEjectSize,trace.size(),utilities.size());
    for(auto it=--(utilities.end()); it!=utilities.begin(); it--) { //TBD last entry?
        assert(it!=utilities.end());
        curEjectSize++;
        if(curEjectSize >= maxEjectSize/2 && (*it != *(--(utilSteps.end())) ) ) {
            utilSteps.push_back(*it);
            //DEBUG
            LOG("utilStep",*it,0,curEjectSize);
            curEjectSize = 0;
        }
    }
    utilSteps.push_back(0); // min util as end
    utilities.clear();
    std::cout << "ejection sets #sets: " << utilSteps.size() << " |set| " << maxEjectSize << "\n";
        

    long double solval;
    size_t effectiveEjectSize;
    std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, int> > lastSeen;
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curNode;
    SmartDigraph::Node prevNode;
    
    // TBD: describe
    long double nonFlexSize = 0;
    std::map<size_t, long double> endOfIntervalSize;

    // LNS iteration steps
    for(size_t k=0; k+2<utilSteps.size(); k++) {
        // set step's util boundaries
        const long double minUtil = utilSteps[k+2];
        const long double maxUtil = utilSteps[k];
        std::cout << "iteration utility min " << minUtil << " max " << maxUtil << "\n";

        // create a graph with just arcs with utility between minUtil and maxUtil
        curNode = g.addNode(); // initial node
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

            // create arcs if in ejection set
            if(isInEjectSet(minUtil, maxUtil, curEntry) ) {
                // second: if there is another request for this object
                if(curEntry.hasNext) {
                    // save prev node as anchor for future arcs
                    prevNode = curNode;
                    lastSeen[std::make_pair(curEntry.id,curEntry.size)]=std::make_pair(i,g.id(prevNode));
                    // create another node, "inner" capacity arc
                    curNode = g.addNode(); // next node
                    curArc = g.addArc(prevNode,curNode);
                    cap[curArc] = cacheSize - std::floor(nonFlexSize); 
                    // DEBUG
                    LOG("iA",i,cap[curArc],0);
                    cost[curArc] = 0;
                }


                //not in ejection set and dvar > 0 -> need to subtract dvar*size for interval's duration
            } else if (curEntry.dvar > 0) { 
                const long double curEffectiveSize = curEntry.size*curEntry.dvar;
                // DEBUG
                LOG("-CS",i,curEffectiveSize,curEntry.dvar);
                assert(cacheSize >= curEffectiveSize);
                nonFlexSize += curEffectiveSize;
                // assert valid flexsize
                endOfIntervalSize.emplace(curEntry.nextSeen,curEffectiveSize);
            }


            // clear all nonFlexSize which are 
            while(endOfIntervalSize.size() > 0 && endOfIntervalSize.begin()->first <= i+1) {
                // DEBUG
                LOG("+CS",i,endOfIntervalSize.begin()->first,endOfIntervalSize.begin()->second);
                nonFlexSize -= endOfIntervalSize.begin()->second;
                endOfIntervalSize.erase(endOfIntervalSize.begin());
            }
        }

        lastSeen.clear();
        endOfIntervalSize.clear();

        // too many nodes

        // solve this MCF
        SmartDigraph::ArcMap<uint64_t> flow(g);
        solval = solveMCF(g, cap, cost, supplies, flow, solverPar);
        std::cout << "sol step cost " << solval << " effES " << effectiveEjectSize <<"\n";

        // write DVAR to trace
        solval = 0;
        for(uint64_t i=0; i<trace.size(); i++) {
            if(trace[i].active) {
                trace[i].dvar = 1.0L - flow[g.arcFromId(trace[i].arcId)]/static_cast<long double>(trace[i].size);
            }
            LOG("dv",i,trace[i].dvar,trace[i].size);
            assert(trace[i].dvar >= 0 && trace[i].dvar<=1);
            solval += trace[i].dvar;
        }

        std::cout << " reqc " << solval << "\n";
    }
    
    return 0;
}



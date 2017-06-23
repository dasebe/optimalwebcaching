#include <fstream>
#include <lemon/lgf_writer.h>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << argv[0] << " traceFile cacheSize maxEjectSize" << std::endl;
        return 1;
    }

    std::cerr << "I'm started and still happy" << std::endl;

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

    // LNS constants
    int64_t totalCapacity = 0;
    int64_t maxFlowAmount = 0;
    double highestCost = -1;

    // dual-feasible solution: all alpha and pi are 0
    uint64_t nodes=0, vertices=0;
    double dualSol = 0;
    SmartDigraph::NodeMap<double> pi(g);
    SmartDigraph::ArcMap<double> alpha(g);

    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodes;
        pi[n] = 0;
        dualSol += pi[n] * supplies[n];
        const int64_t curSuppl = supplies[n];
        if(curSuppl>0) {
            maxFlowAmount += curSuppl;
        }
    }
    std::cerr << nodes << " nodes ";
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++vertices;
        alpha[a] = 0;
        dualSol += alpha[a] * cap[a];
        totalCapacity+=cap[a];
        const double curCost = cost[a];
        if(highestCost > curCost || highestCost==-1) {
            highestCost = curCost;
        }
    }
    std::cerr << vertices << " arcs with dual cost " << dualSol << std::endl;

    for (SmartDigraph::OutArcIt a(g, lastNode); a!=INVALID; ++a) ++nodes;

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc - totalUniqC) {
        maxEjectSize = totalReqc - totalUniqC;
    }
    

    // LNS iteration steps
    // std::unordered_set<int64_t> outerArcs;
    // std::unordered_set<int64_t> innerArcs;
    // std::unordered_map<int64_t,std::vector<int64_t> > lArcs;
    // std::unordered_map<int64_t,std::vector<int64_t> > rArcs;
    std::unordered_set<int64_t> ejectNodes;

    std::unordered_map<int64_t, int64_t> lastSeen;
    // tmp variables
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curNode;
    SmartDigraph::Node prevNode;
    SmartDigraph::Node extraNode;
    int64_t arcId;
    double extraArcCost;
    bool extraArcFlag;


    for(size_t kmin=0; kmin<trace.size(); kmin+=maxEjectSize/2) {
        const size_t kmax = kmin+maxEjectSize;

        // LNS graph structure
        SmartDigraph lnsG; // mcf graph
        extraNode = lnsG.addNode();
        curNode = lnsG.addNode();
        SmartDigraph::ArcMap<int64_t> lnsCap(g); // mcf capacities
        SmartDigraph::ArcMap<double> lnsCost(g); // mcf costs
        SmartDigraph::NodeMap<int64_t> lnsSupplies(g); // mcf demands/supplies

        // get all nodes in the ejection set
        for(uint64_t i=kmin; i<kmax; i++) {
            trEntry & curEntry = trace[i];
            arcId = curEntry.innerArcId;
            const int64_t ejNodeId = g.id(g.source(g.arcFromId(arcId)));
            ejectNodes.insert(ejNodeId);
        }
        // last node
        const int64_t ejNodeId = g.id(g.target(g.arcFromId(arcId)));
        ejectNodes.insert(ejNodeId);

        for(uint64_t i=kmin; i<kmax; i++) {
            trEntry & curEntry = trace[i];

            // handle outer arc
            arcId = curEntry.outerArcId;

            // reset extra arc flag
            extraArcFlag = false;

            // remember current node for outer arcid if in ejection set
            curArc = g.arcFromId(arcId);
            const int64_t iOutNodeId = g.id(g.target(curArc));
            if(ejectNodes.count(iOutNodeId)>0) {
                lastSeen[iOutNodeId] = lnsG.id(curNode);
            } else {
                // not in ejection set
                // c^tile_ij = cij + aij + pi
                extraArcCost = cost[curArc] + alpha[curArc] + pi[g.target(curArc)];
                extraArcFlag = true;
            }

            // handle inner arc
            arcId = curEntry.innerArcId;

            // if inner arc is going out of ejection set
            curArc = g.arcFromId(arcId);
            const int64_t oOutNodeId = g.id(g.target(curArc));
            if(ejectNodes.count(oOutNodeId)==0) {
                // not in ejection set
                // c^tile_ij = cij + aij + pi
                extraArcCost = cost[curArc] + alpha[curArc] + pi[g.target(curArc)];
                extraArcFlag = true;
            }

            // if extra arc flag, create extra arc
            if(extraArcFlag) {
                curArc = g.addArc(extraNode,curNode);
                lnsCost[curArc] = extraArcCost;
                lnsCap[curArc] = maxFlowAmount;
            }

            // iterate over all incoming arcs and derive cost of additional arc
            extraArcFlag = false;
            const SmartDigraph::Node lNode = g.source(g.arcFromId(arcId));
            extraArcCost = highestCost+1;
            for (SmartDigraph::InArcIt a(g, lNode); a!=INVALID; ++a) {
                // c^tile_ij = cij + aij - pi
                const SmartDigraph::Node incArcNodeId = g.source(a);
                const double curCost = cost[a] + alpha[a] - pi[incArcNodeId];
                if(curCost < extraArcCost) {
                    extraArcCost = curCost;
                }
                extraArcFlag = true;
            }
            if(extraArcFlag) {
                // add additional arc to extraNode
                curArc = g.addArc(curNode,extraNode);
                lnsCost[curArc] = extraArcCost;
                lnsCap[curArc] = maxFlowAmount;
            }


            // create all the original graph arcs

            // outer arc
            arcId = curEntry.innerArcId;
            const uint64_t size = curEntry.size;
            if(lastSeen.count(arcId) > 0) {
                const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[arcId]);
                curArc = g.addArc(lastReq,curNode);
                lnsCap[curArc] = size;
                lnsCost[curArc] = 1/static_cast <double>(size);
                lnsSupplies[lastReq] += size;
                lnsSupplies[curNode] -= size;
            }

            // inner arc
            prevNode = curNode;
            curNode = g.addNode(); // next node
            curArc = g.addArc(prevNode,curNode);
            lnsCap[curArc] = cacheSize; 
            lnsCost[curArc] = 0;

            // do we need to do second node??

        }

        std::cerr << "I'm done and still happy" << std::endl;

        SmartDigraph::ArcMap<uint64_t> flow(g);
        double solval = solveMCF(lnsG, lnsCap, lnsCost, lnsSupplies, flow, 4);
        std::cerr << "ExLP" << 4 << " " << cacheSize << " hitc " << totalReqc-totalUniqC-solval << " reqc " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << std::endl;

        return 0;

    }


    return 0;
}

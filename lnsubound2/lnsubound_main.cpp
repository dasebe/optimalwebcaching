#include <fstream>
#include <limits>
#include <lemon/lgf_writer.h>
#include <lemon/core.h>
#include <lemon/bellman_ford.h>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <tuple>
#include <cmath>
#include <chrono>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

#define ZEROEPSILONSET 1e-9
#define ZEROEPSILONCHECK 1e-6

int main(int argc, char* argv[]) {

    if (argc != 5) {
        std::cerr << argv[0] << " traceFile cacheSize maxEjectSize resultPath" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    int64_t cacheSize(atoll(argv[2]));
    uint64_t maxEjectSize(std::stoull(argv[3]));
    std::string resultPath(argv[4]);

    // parse trace file
    std::vector<trEntry> trace;
    OLOG("opening trace",0,0,0);
    uint64_t totalUniqC = parseTraceFile(trace, path);
    uint64_t totalReqc = trace.size();
    OLOG("scanned trace",totalReqc,totalUniqC,0);

    // create mcf instance
    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<long double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supply(g); // mcf demands/supply
    createMCF(g, trace, cacheSize, cap, cost, supply);
    OLOG("finished g",totalReqc,totalUniqC,0);
    
    // LNS constants
    long double highestCost = -1;

    // dual-feasible solution: all alpha and pi are 0
    uint64_t nodeCount = 0, arcCount = 0;
    long double dualSol = 0;
    SmartDigraph::NodeMap<long double> pi(g);
    SmartDigraph::ArcMap<long double> alpha(g);

    SmartDigraph::NodeIt n(g);
    do {
        ++nodeCount;
        pi[n] = 0;
        dualSol += pi[n] * supply[n];
        GLOG("g node",g.id(n),supply[n],0);
    } while (++n != INVALID);
    SmartDigraph::ArcIt a(g);
    do {
        ++arcCount;
        alpha[a] = 0;
        dualSol -= alpha[a] * cap[a];
        const long double curCost = cost[a];
        if(highestCost > curCost || highestCost==-1) {
            highestCost = curCost;
        }
        GLOG("g arc",g.id(g.source(a)),g.id(g.target(a)),cost[a]);
    } while (++a != INVALID);
    OLOG("finished dual init",nodeCount,arcCount,dualSol);

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc) {
        maxEjectSize = totalReqc + 1;
    }
    
    // tmp variables
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curLnsNode;
    SmartDigraph::Node curGNode;
    SmartDigraph::Node extraNode;
    int64_t nodeId;
    int64_t arcId;
    long double extraArcCost;
    bool extraArcFlag;
    uint64_t extraArcCount;
    size_t traceIndex=0, traceHalfIndex = 0;
    long double localDualValue;
    std::chrono::high_resolution_clock::time_point ts, tg, tmcf;

    // iterate over graph again
    for(size_t kmin=trace.size(); kmin>0; kmin=traceHalfIndex) {
        //        OLOG("start config",kmin,traceIndex,maxEjectSize);
        ts = std::chrono::high_resolution_clock::now();

        // LNS graph structure
        SmartDigraph lnsG; // mcf graph
        extraNode = lnsG.addNode();
        nodeCount = 0;
        arcCount = 0;
        extraArcCount = 0;
        // LNS iteration step state
        std::unordered_set<int64_t> ejectNodes;
        // graph min cost flow information
        SmartDigraph::ArcMap<int64_t> lnsCap(lnsG); // mcf capacities
        SmartDigraph::ArcMap<long double> lnsCost(lnsG); // mcf costs
        SmartDigraph::NodeMap<int64_t> lnsSupply(lnsG); // mcf demands/supply
        SmartDigraph::ArcMap<int64_t> lnsFlow(lnsG); // in units of flow
        // temporary pi -- actually unused
        SmartDigraph::NodeMap<long double> lnsPi(lnsG);
        // to map back results to g graph
        SmartDigraph::NodeMap<int64_t> lnsToGNodeId(lnsG);
        SmartDigraph::NodeMap<int64_t> lnsToTraceI(lnsG);
        SmartDigraph::NodeMap<int64_t> gtoLnsNodeId(g);
        lnsToGNodeId[extraNode] = -1;
        SmartDigraph::ArcMap<int64_t> lnsToGArcId(lnsG);

        // mark nodes for the ejection set
        traceIndex = kmin;
        arcId = -1;
        int64_t totalSupply = 0;
        traceHalfIndex = 0;

        while(ejectNodes.size()<=maxEjectSize && traceIndex>0) {
            const trEntry & curEntry = trace[traceIndex-1];
            GLOG("lnsG",traceIndex,curEntry.hasNext,trace.size());
            arcId = curEntry.innerArcId;
            curGNode = INVALID;
            if(curEntry.hasNext) {
                curGNode = g.source(g.arcFromId(arcId));
            } else if (traceIndex==trace.size()) {
                // at end of trace we add the last node
                curGNode = g.target(g.arcFromId(arcId));
            }
            // add node to lns graph
            if(curGNode != INVALID) {
                nodeId = g.id(curGNode);
                ejectNodes.insert(nodeId);
                curLnsNode = lnsG.addNode();
                nodeCount++;
                // copy supply
                lnsSupply[curLnsNode] = supply[curGNode];
                totalSupply += lnsSupply[curLnsNode];
                // establish mapping
                lnsToGNodeId[curLnsNode] = nodeId;
                lnsToTraceI[curLnsNode] = traceIndex;
                gtoLnsNodeId[curGNode] = lnsG.id(curLnsNode);
            }
            // increment index
            traceIndex--;
            // save half time to move forward for loop
            if(ejectNodes.size()==maxEjectSize) {
                traceHalfIndex = traceIndex;
            }
        }
        traceHalfIndex = kmin-maxEjectSize;
        lnsSupply[extraNode] = -totalSupply;
        GLOG("ejectSize",kmin,traceIndex,ejectNodes.size());

        GLOG("trIndex",kmin,traceIndex,0);

        for(int64_t curNodeId: ejectNodes) {
            curGNode = g.nodeFromId(curNodeId);
            curLnsNode = lnsG.nodeFromId(gtoLnsNodeId[curGNode]);

            // INCOMING ARCS
            extraArcFlag = false;
            extraArcCost = highestCost+1;
            SmartDigraph::InArcIt ia(g, curGNode);
            while(ia!=INVALID) {
                const SmartDigraph::Node incSrcNode = g.source(ia);
                const int64_t incSrcNodeId = g.id(incSrcNode);
                GLOG("inc",curNodeId,incSrcNodeId,0);
                if(ejectNodes.count(incSrcNodeId)!=0) {
                    // incoming from ejection set
                    GLOG("add arc ii",incSrcNodeId,curNodeId,cost[ia]);
                    curArc = lnsG.addArc(
                                         lnsG.nodeFromId(gtoLnsNodeId[incSrcNode]),
                                         curLnsNode
                                         );
                    arcCount++;
                    lnsCap[curArc] = cap[ia];
                    lnsCost[curArc] = cost[ia];
                    lnsToGArcId[curArc] = g.id(ia);
                } else {
                    // incoming node not in ejection set
                    // c^tile_ij = cij + aij - pi
                    const long double curCost = cost[ia] + alpha[ia] - pi[incSrcNode];
                    // min over all extra arc costs
                    if(curCost < extraArcCost) {
                        extraArcCost = curCost;
                    }
                    extraArcFlag = true;
                }
                // move loop forward
                ++ia;
            }
            // add incoming arc from EXTRANODE
            if(extraArcFlag) {
                GLOG("add arc ie",lnsToGNodeId[extraNode],curNodeId,extraArcCost);
                curArc = lnsG.addArc(extraNode,curLnsNode);
                extraArcCount++;
                lnsCost[curArc] = extraArcCost;
                lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                lnsToGArcId[curArc] = -1;
            }
        


            // OUTGOING ARCS
            extraArcFlag = false;
            extraArcCost = highestCost+1;
            SmartDigraph::OutArcIt oa(g, curGNode);
            while(oa!=INVALID) {
                const SmartDigraph::Node outTrgtNode = g.target(oa);
                const int64_t outTrgtNodeId = g.id(outTrgtNode);
                GLOG("out",curNodeId,outTrgtNodeId,0);
                if(ejectNodes.count(outTrgtNodeId)!=0) {
                    // constructing both incoming and outgoing arcs would lead to duplicates
                    // outgoing to ejection set
                    // GLOG("add arc oi",curNodeId,outTrgtNodeId,cost[oa]);
                    // curArc = lnsG.addArc(
                    //                      curLnsNode,
                    //                      lnsG.nodeFromId(gtoLnsNodeId[outTrgtNode])
                    //                      );
                    // arcCount++;
                    // lnsCap[curArc] = cap[oa];
                    // lnsCost[curArc] = cost[oa];
                } else {
                    // outgoing node not in ejection set
                    // c^tile_ij = cij + aij + pi
                    GLOG("oe",0,0,0);
                    const long double curExtraArcCost = cost[oa] + alpha[oa] + pi[outTrgtNode];
                    // min over all extra arc costs
                    if(curExtraArcCost < extraArcCost) {
                        extraArcCost = curExtraArcCost;
                    }
                    extraArcFlag = true;
                }
                ++oa;
            }
            // add outgoing arc to EXTRANODE
            if(extraArcFlag) {
                GLOG("add arc oe",curNodeId,lnsToGNodeId[extraNode],extraArcCost);
                curArc = lnsG.addArc(curLnsNode,extraNode);
                extraArcCount++;
                lnsCost[curArc] = extraArcCost;
                lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                lnsToGArcId[curArc] = -1;
            }
        }

        SmartDigraph::ArcIt aIt;

#ifdef GDEBUG
        SmartDigraph::NodeIt nIt(lnsG);
        do {
            GLOG("lns g node",lnsToGNodeId[nIt],lnsSupply[nIt],lnsG.id(nIt));
        } while (++nIt!=INVALID);

        aIt = SmartDigraph::ArcIt(lnsG);
        do {
            GLOG("l arc\t"+std::to_string(lnsToGNodeId[lnsG.source(aIt)]),
                 lnsToGNodeId[lnsG.target(aIt)],
                 lnsCost[aIt],
                 lnsCap[aIt]);
        } while (++aIt!=INVALID);
#endif

        tg = std::chrono::high_resolution_clock::now();
        GLOG("createdLNS",nodeCount,arcCount,extraArcCount);

        // solve the local MCF
        double mcfSol = solveMCF(lnsG, lnsCap, lnsCost, lnsSupply, lnsFlow, 4, lnsPi);
        tmcf = std::chrono::high_resolution_clock::now();

// node potentials based on network simplex

        // PIs are all x(-1) and we need PI[extranode]=0 so subtract it everywhere
        // do this while we compute local dual value (first pi sum)
        localDualValue = 0;
        const long double piOffset = -lnsPi[extraNode];
        SmartDigraph::NodeIt nIt2(lnsG);
        do {
            if(extraNode != nIt2 ) {
                lnsPi[nIt2] = -lnsPi[nIt2] - piOffset;
                localDualValue += lnsPi[nIt2] * lnsSupply[nIt2];
                pi[g.nodeFromId(lnsToGNodeId[nIt2])] = lnsPi[nIt2];
            } else {
                lnsPi[extraNode] = 0;
            }
            PILOG("local pi",lnsToGNodeId[nIt2],lnsPi[nIt2],lnsSupply[nIt2]);
        } while (++nIt2!=INVALID);
        PILOG("local pi sum",localDualValue,0,0);

        aIt = SmartDigraph::ArcIt(lnsG);
        do {
            // calc alpha
            const long double piI = lnsPi[lnsG.source(aIt)];
            const long double piJ = lnsPi[lnsG.target(aIt)];
            const long double cij = lnsCost[aIt];
            long double curAlpha = piI - piJ - cij;
            // correct negative alphas and rounding errors
            if(curAlpha < 0 || std::abs(curAlpha) <= ZEROEPSILONSET) {
                curAlpha = 0.0L;
            }

            PILOG("local alpha\t"+
                  std::to_string(lnsToGNodeId[lnsG.source(aIt)]),
                  lnsToGNodeId[lnsG.target(aIt)],
                  curAlpha,
                  curAlpha * lnsCap[aIt]);
            // update alpha
            if(lnsG.source(aIt)!=extraNode && lnsG.target(aIt)!=extraNode) {
                assert(lnsToGArcId[aIt] != -1);
                alpha[g.arcFromId(lnsToGArcId[aIt])] = curAlpha;
            } else {
                // ensure that alpha for extra arcs are zero
                assert(std::abs(curAlpha) <= ZEROEPSILONCHECK);
            };
            // update local dual objective
            localDualValue -= curAlpha * lnsCap[aIt];
        } while (++aIt!=INVALID);

        //        assert(localDualValue>0);

        // // flow cost setup
        // aIt = SmartDigraph::ArcIt(lnsG);
        // do {
        //         // update dvar/flow mapping in trace vector by calc based on flow in outer arcs
        //         assert(lnsIdToTraceIndex.count(lnsG.id(aIt)) > 0);
        //         trEntry & mapEntry = trace[lnsIdToTraceIndex[lnsG.id(aIt)]];
        //         mapEntry.lcost += static_cast<long double>(lnsFlow[aIt]) * lnsCost[aIt];
        // } while (++aIt!=INVALID);

        OLOG("local dual val\t"+std::to_string(kmin)+"\t"+std::to_string(localDualValue),
             mcfSol,
             std::chrono::duration_cast<std::chrono::duration<long double>>(tg-ts).count(),
             std::chrono::duration_cast<std::chrono::duration<long double>>(tmcf-tg).count()
             );

    }
   

    // long double hitCount = 0;
    // for(auto & it: trace) {
    //     //lcost is not set if no hasNext!!
    //     if(it.hasNext) {
    //         hitCount += 1-it.lcost;
    //     }
    // }



    long double dualVal = 0;
    SmartDigraph::NodeIt nIt(g);
    do {
        dualVal += pi[nIt] * supply[nIt];
        PILOG("global pi",g.id(nIt),pi[nIt],supply[nIt]);
    } while (++nIt!=INVALID);
    PILOG("global pi sum",dualVal,0,0);
    SmartDigraph::ArcIt aIt(g);
    do {
        dualVal -= alpha[aIt] * cap[aIt];
        PILOG("global alpha\t"+
              std::to_string(g.id(g.source(aIt))),
              g.id(g.target(aIt)),
              alpha[aIt],
              alpha[aIt] * cap[aIt]);
    } while (++aIt!=INVALID);
    
    OLOG("global dual val",dualVal,double(totalReqc-totalUniqC)/totalReqc,double(totalReqc-totalUniqC-dualVal)/totalReqc);

    //    std::cerr << "UPB_LNS " << maxEjectSize << " " << cacheSize << " hitc " << hitCount << " reqc " << totalReqc << " OHR " << (static_cast<double>(hitCount))/totalReqc << std::endl;
    //    std::cout << "UPB_LNS " << maxEjectSize << " " << cacheSize << " hitc " << hitCount << " reqc " << totalReqc << " OHR " << (static_cast<double>(hitCount))/totalReqc << std::endl;

    // // output decision variables and utilities
    // std::ofstream resultfile(resultPath);

    // for(auto & it: trace) {
    //     resultfile << it.origTime << " "
    //                << it.id << " " << it.size << " ";
    //     if(it.hasNext) {
    //         resultfile << 1-it.lcost << "\n";
    //     } else {
    //         resultfile << 0 << "\n";
    //     }
    // }


    return 0;
}

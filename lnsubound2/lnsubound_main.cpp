#include <fstream>
#include <limits>
#include <lemon/lgf_writer.h>
#include <lemon/core.h>
#include <lemon/bellman_ford.h>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <queue>
#include <tuple>
#include <cmath>
#include <chrono>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

int main(int argc, char* argv[]) {

    if (argc != 6) {
        std::cerr << argv[0] << " traceFile cacheSize maxEjectSize eps epsAss" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    int64_t cacheSize(std::stoll(argv[2]));
    uint64_t maxEjectSize(std::stoull(argv[3]));
    long double epsilon(std::stold(argv[4]));
    long double epsilonAssert(std::stold(argv[5]));
    //    std::string resultPath(argv[4]);

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

    // ordered list of utilities and check that objects have size less than cache size
    std::multimap<double, size_t> prio;
    for(size_t i=0; i<trace.size(); i++) {
        auto & it = trace[i];
        if(it.size > cacheSize) {
            it.hasNext = false;
        }
        if(it.hasNext) {
            if(it.iLen < 0) {
                OLOG("fail",it.origTime,it.size,it.iLen);
            }
            assert(it.iLen>=0);
            prio.insert(std::make_pair(it.iLen * it.size, i));
        }
    }

    // tmp variables
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curLnsNode;
    SmartDigraph::Node curGNode;
    SmartDigraph::Node extraNode;
    int64_t arcId;
    long double extraArcCost1, extraArcCost2;
    bool extraArcFlag1, extraArcFlag2;
    uint64_t extraArcCount;
    long double localDualValue;
    std::chrono::high_resolution_clock::time_point ts, tg, tmcf;
    std::unordered_set<int64_t> ejectNodes;

    // iterate over graph again
    bool reachedEnd = false;
    auto prioIt=--(prio.end());
    uint64_t processedCount = 0;

    while(!reachedEnd) {
        OLOG("start",prioIt->first,prioIt->second,maxEjectSize);
        ts = std::chrono::high_resolution_clock::now();
        // LNS iteration step state
        ejectNodes.clear();

        // mark nodes for the ejection set
        while(ejectNodes.size() <= maxEjectSize) {
            if(--prioIt == prio.begin()) {
                reachedEnd = true;
                break;
            }
            trEntry & curEntry = trace[prioIt->second];
            //            OLOG("test",i,ejectNodes.size(),0);
            GLOG("lnsG",prioIt->second,curEntry.hasNext,trace.size());
            arcId = curEntry.innerArcId;
            curGNode = INVALID;
            if(curEntry.hasNext) {
                curGNode = g.source(g.arcFromId(arcId));
            } else if (prioIt->second>=trace.size()-1) {
                // at end of trace we add the last node
                curGNode = g.target(g.arcFromId(arcId));
            }
            // add node to lns graph
            if(curGNode != INVALID) {
                const int64_t curId = g.id(curGNode);
                ejectNodes.insert(curId);
                SmartDigraph::InArcIt ia(g, curGNode);
                while(ia!=INVALID) {
                    const int64_t curId = g.id(g.source(ia));
                    ejectNodes.insert(curId);
                    ++ia;
                }
                SmartDigraph::OutArcIt oa(g, curGNode);
                while(oa!=INVALID) {
                    const int64_t curId = g.id(g.target(oa));
                    ejectNodes.insert(curId);
                    ++oa;
                }
            }
            if(!curEntry.processed) {
                curEntry.processed = true;
                processedCount++;
            }
        }
        //        OLOG("ejectSet",kmin,ejectNodes.size(),maxEjectSize);

        // LNS graph structure
        SmartDigraph lnsG; // mcf graph
        extraNode = lnsG.addNode();
        nodeCount = 0;
        arcCount = 0;
        extraArcCount = 0;
        // graph min cost flow information
        SmartDigraph::ArcMap<int64_t> lnsCap(lnsG); // mcf capacities
        SmartDigraph::ArcMap<long double> lnsCost(lnsG); // mcf costs
        SmartDigraph::NodeMap<int64_t> lnsSupply(lnsG); // mcf demands/supply
        SmartDigraph::ArcMap<int64_t> lnsFlow(lnsG); // in units of flow
        // temporary pi -- actually unused
        SmartDigraph::NodeMap<long double> lnsPi(lnsG);
        // to map back results to g graph
        SmartDigraph::NodeMap<int64_t> lnsToGNodeId(lnsG);
        SmartDigraph::NodeMap<int64_t> gtoLnsNodeId(g);
        lnsToGNodeId[extraNode] = -1;
        SmartDigraph::ArcMap<int64_t> lnsToGArcId(lnsG);


        // create lns graph nodes
        int64_t totalSupply = 0;
        for(int64_t nodeId: ejectNodes) {
            curGNode = g.nodeFromId(nodeId);
            curLnsNode = lnsG.addNode();
            nodeCount++;
            // copy supply
            lnsSupply[curLnsNode] = supply[curGNode];
            totalSupply += lnsSupply[curLnsNode];
            // establish mapping
            lnsToGNodeId[curLnsNode] = nodeId;
            gtoLnsNodeId[curGNode] = lnsG.id(curLnsNode);
        }
        lnsSupply[extraNode] = -totalSupply;

        GLOG("ejectSize",processedCount,0,ejectNodes.size());

        for(int64_t curNodeId: ejectNodes) {
            curGNode = g.nodeFromId(curNodeId);
            curLnsNode = lnsG.nodeFromId(gtoLnsNodeId[curGNode]);

            // INCOMING ARCS
            extraArcFlag1 = false;
            extraArcCost1 = highestCost+1;
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
                    GLOG("extraCost ie",curNodeId,incSrcNodeId,cost[ia]);
                    // min over all extra arc costs
                    if(curCost < extraArcCost1) {
                        extraArcCost1 = curCost;
                    }
                    extraArcFlag1 = true;
                }
                // move loop forward
                ++ia;
            }
      


            // OUTGOING ARCS
            extraArcFlag2 = false;
            extraArcCost2 = highestCost+1;
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
                    const long double curExtraArcCost = cost[oa] + alpha[oa] + pi[outTrgtNode];
                    GLOG("extraCost oe",curNodeId,outTrgtNodeId,cost[oa]);
                    // min over all extra arc costs
                    if(curExtraArcCost < extraArcCost2) {
                        extraArcCost2 = curExtraArcCost;
                    }
                    extraArcFlag2 = true;
                }
                ++oa;
            }



            // fix rounding errors in EXTRAarc cost
            if(std::abs(extraArcCost1) <= epsilon) {
                extraArcCost1 = 0;
            }
            if(std::abs(extraArcCost2) <= epsilon) {
                extraArcCost2 = 0;
            }

            if(extraArcCost1<0 && -extraArcCost1 > extraArcCost2) {
                extraArcCost1 = -extraArcCost2;
            }

            if(extraArcCost2<0 && -extraArcCost2 > extraArcCost1) {
                extraArcCost2 = -extraArcCost1;
            }

            if(extraArcFlag1 && extraArcFlag2) {
                if (
                    (extraArcCost1<0 && -extraArcCost1 > extraArcCost2) ||
                    (extraArcCost2<0 && -extraArcCost2 > extraArcCost1)
                    )
                {
                    OLOG("neg cost cycle",curNodeId,extraArcCost1,extraArcCost2);
                }
            }



            // add incoming arc from EXTRANODE
            if(extraArcFlag1) {
                    GLOG("add arc ie",lnsToGNodeId[extraNode],curNodeId,extraArcCost1);
                    curArc = lnsG.addArc(extraNode,curLnsNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost1;
                    if(lnsCost[curArc]<0) {
                        lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                    } else {
                        lnsCap[curArc] = std::numeric_limits<int64_t>::max()-1;
                    }
                    lnsToGArcId[curArc] = -1;
            }




            // add outgoing arc to EXTRANODE
            if(extraArcFlag2) {
                    GLOG("add arc oe",curNodeId,lnsToGNodeId[extraNode],extraArcCost2);
                    curArc = lnsG.addArc(curLnsNode,extraNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost2;
                    if(lnsCost[curArc]<0) {
                        lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                    } else {
                        lnsCap[curArc] = std::numeric_limits<int64_t>::max()-1;
                    }
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
        long double mcfSol = solveMCF(lnsG, lnsCap, lnsCost, lnsSupply, lnsFlow, 4, lnsPi);
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
            if(curAlpha < 0 || std::abs(curAlpha) <= epsilon) {
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
                assert(std::abs(curAlpha) <= epsilonAssert);
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

        // OLOG("local dual val\t"+std::to_string(kmin)+"\t"+std::to_string(localDualValue),
        //      mcfSol,
        //      std::chrono::duration_cast<std::chrono::duration<long double>>(tg-ts).count(),
        //      std::chrono::duration_cast<std::chrono::duration<long double>>(tmcf-tg).count()
        //      );
        OLOG("local dual val\t"+std::to_string(processedCount)+"\t"+std::to_string(localDualValue),
             mcfSol,
             ejectNodes.size(),
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
    
    OLOG("global dual val\t"+std::to_string(dualVal),totalReqc-totalUniqC,static_cast<long double>(totalReqc-totalUniqC)/totalReqc,static_cast<long double>(totalReqc-totalUniqC-dualVal)/totalReqc);

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

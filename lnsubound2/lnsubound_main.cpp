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

#define EPSILON 1e-9
#define EPSILONASSERT 1e-6

using namespace lemon;

void solveLns(std::unordered_set<int64_t> & ejectNodes, SmartDigraph & g, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<long double> & cost, SmartDigraph::NodeMap<int64_t> & supply, SmartDigraph::NodeMap<long double> & pi, SmartDigraph::ArcMap<long double> & alpha) {



    // tmp variables
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curLnsNode;
    SmartDigraph::Node curGNode;
    SmartDigraph::Node extraNode;
    long double extraArcCost1, extraArcCost2;
    bool extraArcFlag1, extraArcFlag2;
    long double localDualValue;
    std::chrono::high_resolution_clock::time_point ts, tg, tmcf;
    ts = std::chrono::high_resolution_clock::now();



    // LNS graph structure
    SmartDigraph lnsG; // mcf graph
    extraNode = lnsG.addNode();
    uint64_t nodeCount = 0, arcCount = 0, extraArcCount = 0;
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

    GLOG("ejectSize",0,0,ejectNodes.size());

    for(int64_t curNodeId: ejectNodes) {
        curGNode = g.nodeFromId(curNodeId);
        curLnsNode = lnsG.nodeFromId(gtoLnsNodeId[curGNode]);

        // INCOMING ARCS
        extraArcFlag1 = false;
        extraArcCost1 = std::numeric_limits<double>::max();
        SmartDigraph::InArcIt ia(g, curGNode);
        while(ia!=INVALID) {
            const SmartDigraph::Node incSrcNode = g.source(ia);
            const int64_t incSrcNodeId = g.id(incSrcNode);
            //                GLOG("inc",curNodeId,incSrcNodeId,0);
            if(ejectNodes.count(incSrcNodeId)!=0) {
                // incoming from ejection set
                //                    GLOG("add arc ii",incSrcNodeId,curNodeId,cost[ia]);
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
                //                    GLOG("extraCost ie",curNodeId,incSrcNodeId,cost[ia]);
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
        extraArcCost2 = std::numeric_limits<double>::max();
        SmartDigraph::OutArcIt oa(g, curGNode);
        while(oa!=INVALID) {
            const SmartDigraph::Node outTrgtNode = g.target(oa);
            const int64_t outTrgtNodeId = g.id(outTrgtNode);
            //                GLOG("out",curNodeId,outTrgtNodeId,0);
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
                //                    GLOG("extraCost oe",curNodeId,outTrgtNodeId,cost[oa]);
                // min over all extra arc costs
                if(curExtraArcCost < extraArcCost2) {
                    extraArcCost2 = curExtraArcCost;
                }
                extraArcFlag2 = true;
            }
            ++oa;
        }



        // fix rounding errors in EXTRAarc cost
        if(std::abs(extraArcCost1) <= EPSILON) {
            extraArcCost1 = 0;
        }
        if(std::abs(extraArcCost2) <= EPSILON) {
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
            //                    GLOG("add arc ie",lnsToGNodeId[extraNode],curNodeId,extraArcCost1);
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
            //                    GLOG("add arc oe",curNodeId,lnsToGNodeId[extraNode],extraArcCost2);
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

    tg = std::chrono::high_resolution_clock::now();
    GLOG("createdLNS",nodeCount,arcCount,extraArcCount);

    // solve the local MCF
    long double mcfSol = solveMCF(lnsG, lnsCap, lnsCost, lnsSupply, lnsFlow, 4, lnsPi);
    tmcf = std::chrono::high_resolution_clock::now();


#ifdef GDEBUG
    // SmartDigraph::NodeIt nIt(lnsG);
    // do {
    //     GLOG("lns g node",lnsToGNodeId[nIt],lnsSupply[nIt],lnsG.id(nIt));
    // } while (++nIt!=INVALID);
        
        
    aIt = SmartDigraph::ArcIt(lnsG);
    do {
        // GLOG("l arc\t"+
        //      std::to_string(lnsToGNodeId[lnsG.source(aIt)])+
        //      "->"+
        //      std::to_string(lnsToGNodeId[lnsG.target(aIt)]),
        //      lnsCost[aIt],
        //      lnsCap[aIt],
        //      lnsFlow[aIt]);
            
    } while (++aIt!=INVALID);
#endif





    // node potentials based on network simplex

    // PIs are all x(-1) and we need PI[extranode]=0 so subtract it everywhere
    // do this while we compute local dual value (first pi sum)

    localDualValue = 0;
    // just some stats
    uint64_t zeroAlpha = 0, zeroPi=0, closeZeroPi = 0, nonzAlpha=0, nonzPi=0;
    const long double piOffset = -lnsPi[extraNode];
    SmartDigraph::NodeIt nIt2(lnsG);
    do {
        if(extraNode != nIt2 ) {
            lnsPi[nIt2] = -lnsPi[nIt2] - piOffset;
            localDualValue += lnsPi[nIt2] * lnsSupply[nIt2];
            pi[g.nodeFromId(lnsToGNodeId[nIt2])] = lnsPi[nIt2];
            // extra stats
            if(lnsPi[nIt2] == 0)
                zeroPi ++;
            else if(std::abs(lnsPi[nIt2]) <= EPSILON) {
                closeZeroPi ++;
            } else
                nonzPi++;
        } else {
            lnsPi[extraNode] = 0;
        }
        PILOG("local pi",lnsToGNodeId[nIt2],lnsPi[nIt2],lnsSupply[nIt2]);
    } while (++nIt2!=INVALID);
    PILOG("local pi sum",localDualValue,0,0);
    OLOG("PIstats",zeroPi,closeZeroPi,nonzPi);

    aIt = SmartDigraph::ArcIt(lnsG);
    do {
        // calc alpha
        const long double piI = lnsPi[lnsG.source(aIt)];
        const long double piJ = lnsPi[lnsG.target(aIt)];
        const long double cij = lnsCost[aIt];
        long double curAlpha = piI - piJ - cij;
        // correct negative alphas and rounding errors
        if(curAlpha < 0 || std::abs(curAlpha) <= EPSILON) {
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
            if(std::abs(curAlpha) <= EPSILON) {
                zeroAlpha++;
            } else
                nonzAlpha++;
        } else {
            // ensure that alpha for extra arcs are zero
            assert(std::abs(curAlpha) <= EPSILONASSERT);
        };
        // update local dual objective
        localDualValue -= curAlpha * lnsCap[aIt];
    } while (++aIt!=INVALID);

    OLOG("Alphastats",zeroAlpha,0,nonzAlpha);

    OLOG("local dual val\t"+std::to_string(0)+"\t"+std::to_string(localDualValue),
         mcfSol,
         ejectNodes.size(),
         std::chrono::duration_cast<std::chrono::duration<long double>>(tmcf-tg).count()
         );
}



int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << argv[0] << " traceFile cacheSize maxEjectSize" << std::endl;
        return 1;
    }

    std::string path(argv[1]);
    int64_t cacheSize(std::stoll(argv[2]));
    uint64_t maxEjectSize(std::stoull(argv[3]));
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
    SmartDigraph::NodeMap<int64_t> traceIndex(g);
    createMCF(g, trace, cacheSize, cap, cost, supply, traceIndex);
    OLOG("finished g",totalReqc,totalUniqC,0);
    
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
        //GLOG("g node",g.id(n),supply[n],0);
    } while (++n != INVALID);
    SmartDigraph::ArcIt a(g);
    do {
        ++arcCount;
        alpha[a] = 0;
        dualSol -= alpha[a] * cap[a];
        //GLOG("g arc",g.id(g.source(a)),g.id(g.target(a)),cost[a]);
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
            const double curPrio = (it.iLen * it.size); //1.0/
            prio.insert(std::make_pair(curPrio, i));
        }
    }
    

    SmartDigraph::Node curGNode;
    std::unordered_set<int64_t> ejectNodes;

    // iterate over graph again
    bool reachedEnd = false;
    size_t kmin=trace.size()-1;

    while(!reachedEnd) {
        OLOG("start1",kmin,0,maxEjectSize);
        // LNS iteration step state
        ejectNodes.clear();

        while(kmin-- >0 && ejectNodes.size() <= maxEjectSize) {
            if(kmin==0)
                reachedEnd = true;
            const trEntry & curEntry = trace[kmin];
            const int64_t arcId = curEntry.innerArcId;
            curGNode = INVALID;
            if(curEntry.hasNext) {
                curGNode = g.source(g.arcFromId(arcId));
            } else if (kmin==trace.size()) {
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
        }

        solveLns(ejectNodes, g, cap, cost, supply, pi, alpha);

    }

    long double dualVal = 0;
    SmartDigraph::NodeIt nIt(g);
    do {
        dualVal += pi[nIt] * supply[nIt];
    } while (++nIt!=INVALID);
    SmartDigraph::ArcIt aIt(g);
    do {
        dualVal -= alpha[aIt] * cap[aIt];
    } while (++aIt!=INVALID);
    
    OLOG("global dual val1\t"+std::to_string(dualVal),totalReqc-totalUniqC,static_cast<long double>(totalReqc-totalUniqC)/totalReqc,static_cast<long double>(totalReqc-totalUniqC-dualVal)/totalReqc);



    reachedEnd = false;
    auto prioIt=--(prio.end());
    uint64_t processedCount = 0;

    while(!reachedEnd) {
        OLOG("start2",prioIt->first,prioIt->second,processedCount);
        // LNS iteration step state
        ejectNodes.clear();

        // mark nodes for the ejection set
        while(ejectNodes.size() <= maxEjectSize) {
            if(--prioIt == prio.begin()) {
                reachedEnd = true;
                break;
            }
            trEntry & curEntry = trace[prioIt->second];
            //            GLOG("lnsG",prioIt->second,curEntry.hasNext,trace.size());
            const int64_t arcId = curEntry.innerArcId;
            curGNode = INVALID;
            if(curEntry.hasNext) {
                curGNode = g.source(g.arcFromId(arcId));
            } else if (prioIt->second>=trace.size()-1) {
                // at end of trace we add the last node
                curGNode = g.target(g.arcFromId(arcId));
            }
            // add node to lns graph
            if(curGNode != INVALID) {
                std::queue<int64_t> neighborNodes;
                neighborNodes.push(g.id(curGNode));
                for(unsigned int d=0; d<100; d++) {
                    const int64_t curId = neighborNodes.front();
                    curGNode = g.nodeFromId(curId);
                    neighborNodes.pop();
                    ejectNodes.insert(curId);
                    SmartDigraph::InArcIt ia(g, curGNode);
                    while(ia!=INVALID) {
                        const int64_t newId = g.id(g.source(ia));
                        neighborNodes.push(newId);
                        ++ia;
                    }
                    SmartDigraph::OutArcIt oa(g, curGNode);
                    while(oa!=INVALID) {
                        const int64_t newId = g.id(g.target(oa));
                        neighborNodes.push(newId);
                        ++oa;
                    }
                }
            }
            if(!curEntry.processed) {
                curEntry.processed = true;
                processedCount++;
            }
            
        }
        if(ejectNodes.size()==0 && reachedEnd) {
            // we're done
            break;
        }

        solveLns(ejectNodes, g, cap, cost, supply, pi, alpha);


    dualVal = 0;
    nIt = SmartDigraph::NodeIt(g);
    do {
        dualVal += pi[nIt] * supply[nIt];
        PILOG("global pi",g.id(nIt),pi[nIt],supply[nIt]);
    } while (++nIt!=INVALID);
    PILOG("global pi sum",dualVal,0,0);
    aIt = SmartDigraph::ArcIt(g);
    do {
        dualVal -= alpha[aIt] * cap[aIt];
        PILOG("global alpha\t"+
              std::to_string(g.id(g.source(aIt))),
              g.id(g.target(aIt)),
              alpha[aIt],
              alpha[aIt] * cap[aIt]);
    } while (++aIt!=INVALID);
    
    OLOG("global dual val\t"+std::to_string(dualVal),totalReqc-totalUniqC,static_cast<long double>(totalReqc-totalUniqC)/totalReqc,static_cast<long double>(totalReqc-totalUniqC-dualVal)/totalReqc);


    }
   

    // long double hitCount = 0;
    // for(auto & it: trace) {
    //     //lcost is not set if no hasNext!!
    //     if(it.hasNext) {
    //         hitCount += 1-it.lcost;
    //     }
    // }




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

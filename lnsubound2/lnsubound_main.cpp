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

#define ZEROEPSILON 1e-8

typedef BellmanFord<SmartDigraph, SmartDigraph::ArcMap<double>> BellSolveType;

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
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supply(g); // mcf demands/supply
    createMCF(g, trace, cacheSize, cap, cost, supply);
    OLOG("finished g",totalReqc,totalUniqC,0);
    
    // LNS constants
    double highestCost = -1;

    // dual-feasible solution: all alpha and pi are 0
    uint64_t nodeCount = 0, arcCount = 0;
    double dualSol = 0;
    SmartDigraph::NodeMap<double> pi(g);
    SmartDigraph::ArcMap<double> alpha(g);

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
        const double curCost = cost[a];
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
    SmartDigraph::Node curNode;
    SmartDigraph::Node prevNode;
    SmartDigraph::Node extraNode;
    int64_t arcId;
    double extraArcCost;
    bool extraArcFlag;
    size_t traceIndex=0, traceHalfIndex = 0;
    int64_t ejNodeIdT;
    long double localDualValue, globalDualValue=0;
    std::chrono::high_resolution_clock::time_point ts, tg, tmcf, tdone;
    OLOG("bLoop",0,0,0);

    // iterate over graph again
    for(size_t kmin=0; kmin<trace.size(); kmin=traceHalfIndex+1) {
        OLOG("start config",kmin,traceIndex,maxEjectSize);
        ts = std::chrono::high_resolution_clock::now();

        // LNS graph structure
        SmartDigraph lnsG; // mcf graph
        extraNode = lnsG.addNode();
        curNode = lnsG.addNode();
        prevNode = curNode;
        nodeCount = 0;
        arcCount = 0;
        uint64_t extraArcCount = 0;
        // LNS iteration step state
        std::unordered_set<int64_t> ejectNodes;
        std::unordered_map<std::pair<uint64_t, uint64_t>, std::tuple<int64_t, int64_t, size_t> > lastSeen;
        // graph min cost flow information
        SmartDigraph::ArcMap<int64_t> lnsCap(lnsG); // mcf capacities
        SmartDigraph::ArcMap<double> lnsCost(lnsG); // mcf costs
        SmartDigraph::NodeMap<int64_t> lnsSupply(lnsG); // mcf demands/supply
        SmartDigraph::ArcMap<int64_t> lnsFlow(lnsG); // in units of flow
        // temporary pi -- actually unused
        SmartDigraph::NodeMap<double> lnsPi(lnsG);
        // to map back results to g graph
        SmartDigraph::NodeMap<int64_t> lnsToGNodeId(lnsG);
        lnsToGNodeId[extraNode] = -1;
        SmartDigraph::ArcMap<int64_t> lnsToGArcId(lnsG);
        std::unordered_map<int64_t, size_t > lnsIdToTraceIndex;

        // mark nodes for the ejection set
        traceIndex = kmin;
        arcId = -1;
        while(ejectNodes.size()<=maxEjectSize && traceIndex<trace.size()) {
            const trEntry & curEntry = trace[traceIndex];
            GLOG("lnsG",traceIndex,curEntry.hasNext,trace.size());
            
            if(curEntry.hasNext) {
                arcId = curEntry.innerArcId;
                const int64_t ejNodeIdS = g.id(g.source(g.arcFromId(arcId)));
                ejectNodes.insert(ejNodeIdS);
                // if at end of trace, we also must add the last node so save it
                ejNodeIdT = g.id(g.target(g.arcFromId(arcId)));
            }
            // at end of trace, so add last node
            if(traceIndex==trace.size()-1) {
                ejectNodes.insert(ejNodeIdT);
            }
            // calc max flow amount as sum of pos supply
            traceIndex++;
            // save half time to move forward for loop
            if(traceIndex-kmin<=maxEjectSize) {
                traceHalfIndex = traceIndex;
            }
        }
        assert(arcId >= 0);
        GLOG("ejectSize",kmin,traceIndex,ejectNodes.size());

        GLOG("trIndex",kmin,traceIndex,0);

        uint64_t localUniqCount = 0;
        for(size_t i=kmin; i<traceIndex; i++) {
            const trEntry & curEntry = trace[i];
            auto curIdSize = std::make_pair(curEntry.id,curEntry.size);
            //            LLOG("curE",i,g.id(g.source(g.arcFromId(curEntry.outerArcId))),g.id(g.target(g.arcFromId(curEntry.outerArcId))));
            GLOG("entry",i,curEntry.hasNext,trace.size());

            // create outer arc, if one is saved for this curIdSize
            if(lastSeen.count(curIdSize) > 0) { // 
                GLOG("lastSeenNode",lnsToGNodeId[curNode],i,trace.size());
                const uint64_t size = curEntry.size;
                const SmartDigraph::Node lastReq = lnsG.nodeFromId(std::get<0>(lastSeen[curIdSize]));
                curArc = lnsG.addArc(lastReq,curNode);
                arcCount++;
                lnsCap[curArc] = size;
                lnsCost[curArc] = 1/static_cast <double>(size);
                // assert that no link to extranode or nonexisting node here
                assert(std::get<1>(lastSeen[curIdSize]) != -1);
                lnsToGArcId[curArc] = std::get<1>(lastSeen[curIdSize]);
                lnsIdToTraceIndex[lnsG.id(curArc)] = std::get<2>(lastSeen[curIdSize]);
                LLOG("lastseen--",i,curEntry.id,curEntry.size);
                lastSeen.erase(curIdSize);
            }
            else {
                localUniqCount++;
                LLOG("lastseen()",i,curEntry.id,curEntry.size);
            }

            // if this curIdSize has another request in the future, create graph node
            if(curEntry.hasNext) {
                // handle outer arc
                arcId = curEntry.outerArcId;

                // reset extra arc flags
                extraArcFlag = false;
                extraArcCost = highestCost+1;

                // remember current node for outer arcid if in ejection set
                curArc = g.arcFromId(arcId);
                const int64_t iOutNodeId = g.id(g.target(curArc));
                if(ejectNodes.count(iOutNodeId)>0) {
                    GLOG("o node in E",iOutNodeId,0,0);
                    LLOG("lastseen++",i,curEntry.id,curEntry.size);
                    lastSeen[curIdSize] = std::make_tuple(lnsG.id(curNode),g.id(curArc),i);
                } else {
                    GLOG("o not in E",iOutNodeId,0,0);
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
                    GLOG("i not in E",iOutNodeId,0,0);
                    // not in ejection set
                    // c^tile_ij = cij + aij + pi
                    const double curExtraArcCost = cost[curArc] + alpha[curArc] + pi[g.target(curArc)];
                    // min over all extra arc costs
                    if(curExtraArcCost < extraArcCost) {
                        extraArcCost = curExtraArcCost;
                    }
                    extraArcFlag = true;
                } else {
                    GLOG("i node in E",iOutNodeId,0,0);
                }

                // if extra arc flag, create extra arc
                if(extraArcFlag) {
                    GLOG("cur->extra",lnsG.id(curNode),lnsG.id(extraNode),0);
                    curArc = lnsG.addArc(curNode,extraNode); //curArc = lnsG.addArc(extraNode,curNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost;
                    lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                    lnsIdToTraceIndex[lnsG.id(curArc)] = i;
                }

                // add inner node
                prevNode = curNode;
                curNode = lnsG.addNode(); // next node
                nodeCount++;

                // copy over supply from g
                lnsSupply[prevNode] = supply[g.source(g.arcFromId(arcId))];
                lnsSupply[curNode] = supply[g.target(g.arcFromId(arcId))];

                // save mapping of curNode to g graph node
                lnsToGNodeId[prevNode] = g.id(g.source(g.arcFromId(arcId)));
                lnsToGNodeId[curNode] = g.id(g.target(g.arcFromId(arcId)));

                GLOG("added node",lnsToGNodeId[curNode],0,lnsToGNodeId[prevNode]);

                // create inner arc
                curArc = lnsG.addArc(prevNode,curNode);
                arcCount++;
                lnsCap[curArc] = cacheSize; 
                lnsCost[curArc] = 0;
                assert(curEntry.innerArcId!=-1);
                lnsToGArcId[curArc] = curEntry.innerArcId;
                lnsIdToTraceIndex[lnsG.id(curArc)] = i;
                LLOG("curEA",i,curEntry.id,lnsG.id(curNode));
            }

            // iterate over all incoming arcs of current node
            // at end of trace, so add last node
            if(curEntry.hasNext) {
                extraArcFlag = false;
                extraArcCost = highestCost+1;
                arcId = curEntry.innerArcId;
                const SmartDigraph::Node lNode = g.source(g.arcFromId(arcId));
                SmartDigraph::InArcIt a(g, lNode);
                while(a!=INVALID) {
                    const SmartDigraph::Node incArcNode = g.source(a);
                    const int64_t incArcNodeId = g.id(incArcNode);
                    GLOG("incArcs iter",i,g.id(lNode),incArcNodeId);
                    if(ejectNodes.count(incArcNodeId)==0) {
                        // not in ejection set
                        // c^tile_ij = cij + aij - pi
                        const double curCost = cost[a] + alpha[a] - pi[incArcNode];
                        if(curCost < extraArcCost) {
                            extraArcCost = curCost;
                        }
                        extraArcFlag = true;
                        GLOG("incArcs ext",i,g.id(lNode),incArcNodeId);
                    }
                    // move loop forward
                    ++a;
                }
                OLOG("done",0,0,0);
                if(extraArcFlag) {
                    //                const SmartDigraph::Node targetNode = 
                    GLOG("extra->cur",i,lnsToGNodeId[extraNode],lnsToGNodeId[prevNode]);
                    // add additional arc to extraNode
                    curArc = lnsG.addArc(extraNode,prevNode); //                     curArc = lnsG.addArc(curNode,extraNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost;
                    lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                    lnsIdToTraceIndex[lnsG.id(curArc)] = i;
                }
            }

            if(i==trace.size()-1) {
                OLOG("last one check",0,0,0);
                extraArcFlag = false;
                extraArcCost = highestCost+1;
                arcId = curEntry.innerArcId;
                const SmartDigraph::Node lNode = g.target(g.arcFromId(arcId));
                SmartDigraph::InArcIt a(g, lNode);
                do {
                    const SmartDigraph::Node incArcNode = g.source(a);
                    const int64_t incArcNodeId = g.id(incArcNode);
                    GLOG("incArcs iter",i,g.id(lNode),incArcNodeId);
                    if(ejectNodes.count(incArcNodeId)==0) {
                        // not in ejection set
                        // c^tile_ij = cij + aij - pi
                        const double curCost = cost[a] + alpha[a] - pi[incArcNode];
                        if(curCost < extraArcCost) {
                            extraArcCost = curCost;
                        }
                        extraArcFlag = true;
                        GLOG("incArcs ext",i,g.id(lNode),incArcNodeId);
                    }
                }  while(++a!=INVALID);
                if(extraArcFlag) {
                    //                const SmartDigraph::Node targetNode = 
                    GLOG("extra->cur",i,lnsToGNodeId[extraNode],lnsToGNodeId[curNode]);
                    // add additional arc to extraNode
                    curArc = lnsG.addArc(extraNode,curNode); //                     curArc = lnsG.addArc(curNode,extraNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost;
                    lnsCap[curArc] = std::numeric_limits<int64_t>::max();
                    lnsIdToTraceIndex[lnsG.id(curArc)] = i;
                }
            }



        }

        // check supply
        int64_t totalSupply = 0;
        SmartDigraph::NodeIt nIt(lnsG);
        do {
            if(extraNode != nIt ) {
                totalSupply += lnsSupply[nIt];
            }
        } while (++nIt!=INVALID);
        lnsSupply[extraNode] = -totalSupply;
        OLOG("totalSupply",totalSupply,lnsSupply[extraNode],0);

#ifdef LDEBUG
        for(auto it: lastSeen) {
            LLOG("left-over lastSeen",it.first.first,it.first.second,std::get<0>(it.second));
        }
#endif
        assert(lastSeen.size()==0);

        GLOG("extranode lastnode",lnsG.id(extraNode),lnsG.id(curNode),0);
        SmartDigraph::ArcIt aIt;

#ifdef GDEBUG
        nIt = SmartDigraph::NodeIt(lnsG);
        do {
            GLOG("lns g node",lnsToGNodeId[nIt],lnsSupply[nIt],0);
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
        double mcfsol = solveMCF(lnsG, lnsCap, lnsCost, lnsSupply, lnsFlow, 4, lnsPi);
        OLOG("MCF sol val",mcfsol,0,0);
        tmcf = std::chrono::high_resolution_clock::now();

// node potentials based on network simplex
        // PIs are all x(-1) and we need PI[extranode]=0 so subtract it everywhere
        // do this while we compute local dual value (first pi sum)
        localDualValue = 0;
        const double piOffset = -lnsPi[extraNode];
        SmartDigraph::NodeIt nIt2(lnsG);
        do {
            if(extraNode != nIt2 ) {
                lnsPi[nIt2] = -lnsPi[nIt2] - piOffset;
                localDualValue += lnsPi[nIt2] * lnsSupply[nIt2];
                pi[g.nodeFromId(lnsToGNodeId[nIt2])] = lnsPi[nIt2];
            } else {
                lnsPi[extraNode] = 0;
            }
            PILOG("pi old val",lnsToGNodeId[nIt2],lnsPi[nIt2],lnsSupply[nIt2]);
        } while (++nIt2!=INVALID);
        OLOG("OLD pi sum",localDualValue,0,0);

        
        SmartDigraph::ArcMap<double> lnsAlpha(lnsG);
        aIt = SmartDigraph::ArcIt(lnsG);
        do {
            const double piI = lnsPi[lnsG.source(aIt)];
            const double piJ = lnsPi[lnsG.target(aIt)];
            const double cij = lnsCost[aIt];
            const double curAlpha = piI - piJ - cij;
            // calculate alpha
            if(curAlpha > 0)
                lnsAlpha[aIt] = curAlpha;
            else
                lnsAlpha[aIt] = 0;
            PILOG("alpha old val",lnsToGNodeId[lnsG.source(aIt)],lnsToGNodeId[lnsG.target(aIt)],lnsAlpha[aIt]); //,lnsCap[aIt],lnsAlpha[aIt] * lnsCap[aIt]);
            // update alpha
            if(lnsG.source(aIt)!=extraNode && lnsG.target(aIt)!=extraNode) {
                const int64_t sId = lnsToGNodeId[lnsG.source(aIt)], tId = lnsToGNodeId[lnsG.target(aIt)];
                alpha[findArc<SmartDigraph>(g,g.nodeFromId(sId),g.nodeFromId(tId))] = lnsAlpha[aIt];
            } else {
                // ensure that alpha for extra arcs are zero
                if(std::abs(lnsAlpha[aIt]) <= ZEROEPSILON) {
                    lnsAlpha[aIt] = 0;
                } else {
                    assert(lnsAlpha[aIt] == 0);
                }
            };
            // update local dual objective
            localDualValue -= lnsAlpha[aIt] * lnsCap[aIt];
        } while (++aIt!=INVALID);

        //        assert(localDualValue>0);

        // flow cost setup
        aIt = SmartDigraph::ArcIt(lnsG);
        do {
                // update dvar/flow mapping in trace vector by calc based on flow in outer arcs
                assert(lnsIdToTraceIndex.count(lnsG.id(aIt)) > 0);
                trEntry & mapEntry = trace[lnsIdToTraceIndex[lnsG.id(aIt)]];
                mapEntry.lcost += static_cast<double>(lnsFlow[aIt]) * lnsCost[aIt];
        } while (++aIt!=INVALID);

        globalDualValue += localDualValue;

        tdone = std::chrono::high_resolution_clock::now();
        OLOG("NEW dual val",localDualValue,0,localUniqCount);
        OLOG("timings",
             std::chrono::duration_cast<std::chrono::duration<double>>(tg-ts).count(),
             std::chrono::duration_cast<std::chrono::duration<double>>(tmcf-tg).count(),
             std::chrono::duration_cast<std::chrono::duration<double>>(tdone-tmcf).count()
             );

    }

    double hitCount = 0;
    for(auto & it: trace) {
        //lcost is not set if no hasNext!!
        if(it.hasNext) {
            hitCount += 1-it.lcost;
        }
    }

    OLOG("final dual Val",globalDualValue,hitCount,totalReqc-totalUniqC-globalDualValue);

    std::cerr << "UPB_LNS " << maxEjectSize << " " << cacheSize << " hitc " << hitCount << " reqc " << totalReqc << " OHR " << (static_cast<double>(hitCount))/totalReqc << std::endl;
    std::cout << "UPB_LNS " << maxEjectSize << " " << cacheSize << " hitc " << hitCount << " reqc " << totalReqc << " OHR " << (static_cast<double>(hitCount))/totalReqc << std::endl;

    // output decision variables and utilities
    std::ofstream resultfile(resultPath);

    for(auto & it: trace) {
        resultfile << it.origTime << " "
                   << it.id << " " << it.size << " ";
        if(it.hasNext) {
            resultfile << 1-it.lcost << "\n";
        } else {
            resultfile << 0 << "\n";
        }
    }


    return 0;
}

#include <fstream>
#include <lemon/lgf_writer.h>
#include <lemon/bellman_ford.h>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <math.h>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

using namespace lemon;

typedef BellmanFord<SmartDigraph, SmartDigraph::ArcMap<double>> BellSolveType;

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << argv[0] << " traceFile cacheSize maxEjectSize" << std::endl;
        return 1;
    }

    std::cerr << "I'm started and still happy" << std::endl;

    std::string path(argv[1]);
    int64_t cacheSize(atoll(argv[2]));
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
    createMCF(g, trace, cacheSize, cap, cost, supplies);
    
    // LNS constants
    int64_t totalCapacity = 0;
    int64_t maxFlowAmount = 0;
    double highestCost = -1;

    // dual-feasible solution: all alpha and pi are 0
    uint64_t nodeCount = 0, arcCount = 0;
    double dualSol = 0;
    SmartDigraph::NodeMap<double> pi(g);
    SmartDigraph::ArcMap<double> alpha(g);

    
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodeCount;
        pi[n] = 0;
        dualSol += pi[n] * supplies[n];
        //        LOG("g node",g.id(n),0,0);
    }
    SmartDigraph::ArcIt a(g);
    do {
        ++arcCount;
        alpha[a] = 0;
        dualSol -= alpha[a] * cap[a];
        totalCapacity+=cap[a];
        const double curCost = cost[a];
        if(highestCost > curCost || highestCost==-1) {
            highestCost = curCost;
        }
        //        LOG("g arc",g.id(g.source(a)),g.id(g.target(a)),cost[a]);
    } while (++a != INVALID);
    std::cerr << "created graph with "
              << nodeCount << " nodes "
              << arcCount << " arcs with dual cost "
              << dualSol << std::endl << std::endl;

    // max ejection size mustn't be larger than actual trace
    if(maxEjectSize > totalReqc - totalUniqC) {
        maxEjectSize = totalReqc - totalUniqC + 1;
    }
    
    // tmp variables
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curNode;
    SmartDigraph::Node prevNode;
    SmartDigraph::Node extraNode;
    int64_t arcId;
    double extraArcCost;
    bool extraArcFlag;
    size_t traceIndex;

    // todo fixme: kmin
    for(size_t kmin=0; kmin<trace.size(); kmin=traceIndex - maxEjectSize/2) {
        // LNS graph structure
        SmartDigraph lnsG; // mcf graph
        extraNode = lnsG.addNode();
        curNode = lnsG.addNode();
        nodeCount = 0;
        arcCount = 0;
        uint64_t extraArcCount = 0;
        // LNS iteration step state
        std::unordered_set<int64_t> ejectNodes;
        std::unordered_map<std::pair<uint64_t, uint64_t>, int64_t > lastSeen;
        //std::unordered_map<int64_t, int64_t> lastSeen;
        // graph min cost flow information
        SmartDigraph::ArcMap<int64_t> lnsCap(lnsG); // mcf capacities
        SmartDigraph::ArcMap<double> lnsCost(lnsG); // mcf costs
        SmartDigraph::NodeMap<int64_t> lnsSupplies(lnsG); // mcf demands/supplies
        SmartDigraph::ArcMap<int64_t> lnsFlow(lnsG); // in units of flow
        // temporary pi -- actually unused
        SmartDigraph::NodeMap<double> lnsPi(lnsG);
        // to map back results to g graph
        SmartDigraph::NodeMap<int64_t> lnsToGNodeId(lnsG);
        lnsToGNodeId[extraNode] = -1;
        SmartDigraph::ArcMap<int64_t> lnsToGArcId(lnsG);

        // mark nodes for the ejection set
        traceIndex = kmin;
        arcId = -1;
        maxFlowAmount = 0;
        while(ejectNodes.size()<=maxEjectSize && traceIndex<trace.size()) {
            const trEntry & curEntry = trace[traceIndex];
            //            LOG("lnsG",traceIndex,curEntry.hasNext,trace.size());
            if(curEntry.hasNext) {
                arcId = curEntry.innerArcId;
                const int64_t ejNodeIdS = g.id(g.source(g.arcFromId(arcId)));
                ejectNodes.insert(ejNodeIdS);
                const int64_t ejNodeIdT = g.id(g.target(g.arcFromId(arcId)));
                ejectNodes.insert(ejNodeIdT);
                // calc max flow amount as sum of pos supplies
                maxFlowAmount += curEntry.size;
            }
            traceIndex++;
        }
        assert(maxFlowAmount>0);
        assert(arcId >= 0);
        LOG("ejS",kmin,maxEjectSize,ejectNodes.size());

        /*        for(auto it: ejectNodes) {
            LOG("lns g node",it,0,0);
            }*/

        //        LOG("trIndex",kmin,traceIndex,0);

        for(size_t i=kmin; i<traceIndex; i++) {
            const trEntry & curEntry = trace[i];
            auto curIdSize = std::make_pair(curEntry.id,curEntry.size);
            //            LOG("curE",i,curEntry.id,curEntry.hasNext);

            // create outer arc, if one is saved for this curIdSize
            if(lastSeen.count(curIdSize) > 0) {
                const uint64_t size = curEntry.size;
                const SmartDigraph::Node lastReq = lnsG.nodeFromId(lastSeen[curIdSize]);
                //                LOG("curEL",i,lnsG.id(lastReq),lnsG.id(curNode));
                curArc = lnsG.addArc(lastReq,curNode);
                arcCount++;
                lnsCap[curArc] = size;
                lnsCost[curArc] = 1/static_cast <double>(size);
                lnsSupplies[lastReq] += size;
                lnsSupplies[curNode] -= size;
                lnsToGArcId[curArc] = curEntry.outerArcId;
                //                LOG("lastseen--",i,lnsToGNodeId[lastReq],lnsToGNodeId[curNode]);
                lastSeen.erase(curIdSize);
            } /*else {
                LOG("lastseen()",i,curEntry.id,curEntry.size);
                }*/
            
            // if this curIdSize has another request in the future, create graph node
            if(curEntry.hasNext) {
                // reset extra arc flag
                extraArcFlag = false;
                extraArcCost = highestCost+1;

                // handle outer arc
                arcId = curEntry.outerArcId;

                // remember current node for outer arcid if in ejection set
                curArc = g.arcFromId(arcId);
                const int64_t iOutNodeId = g.id(g.target(curArc));
                if(ejectNodes.count(iOutNodeId)>0) {
                    //                    LOG("lastseen++",i,lnsToGNodeId[curNode],curEntry.id);
                    lastSeen[curIdSize] = lnsG.id(curNode);
                } else {
                    LOG("not in ejection set",iOutNodeId,0,0);
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
                    LOG("not in ejection set",oOutNodeId,0,0);
                    // not in ejection set
                    // c^tile_ij = cij + aij + pi
                    const double curExtraArcCost = cost[curArc] + alpha[curArc] + pi[g.target(curArc)];
                    // min over all extra arc costs
                    if(curExtraArcCost < extraArcCost) {
                        extraArcCost = curExtraArcCost;
                    }
                    extraArcFlag = true;
                }

                // if extra arc flag, create extra arc
                if(extraArcFlag) {
                    curArc = lnsG.addArc(extraNode,curNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost;
                    lnsCap[curArc] = maxFlowAmount;
                    // reset extra arc flag
                    extraArcFlag = false;
                    extraArcCost = highestCost+1;
                }

                // iterate over all incoming arcs and derive cost of additional arc
                const SmartDigraph::Node lNode = g.source(g.arcFromId(arcId));
                for (SmartDigraph::InArcIt a(g, lNode); a!=INVALID; ++a) {
                    const SmartDigraph::Node incArcNode = g.source(a);
                    const int64_t incArcNodeId = g.id(incArcNode);
                    if(ejectNodes.count(incArcNodeId)==0) {
                        // not in ejection set
                        // c^tile_ij = cij + aij - pi
                        const double curCost = cost[a] + alpha[a] - pi[incArcNode];
                        if(curCost < extraArcCost) {
                            extraArcCost = curCost;
                        }
                        extraArcFlag = true;
                    }
                }
                if(extraArcFlag) {
                    // add additional arc to extraNode
                    curArc = lnsG.addArc(curNode,extraNode);
                    extraArcCount++;
                    lnsCost[curArc] = extraArcCost;
                    lnsCap[curArc] = maxFlowAmount;
                }

                // add inner node
                prevNode = curNode;
                curNode = lnsG.addNode(); // next node
                nodeCount++;

                // save mapping of curNode to g graph node
                lnsToGNodeId[prevNode] = g.id(g.source(g.arcFromId(arcId)));
                lnsToGNodeId[curNode] = g.id(g.target(g.arcFromId(arcId)));

                // create inner arc
                curArc = lnsG.addArc(prevNode,curNode);
                arcCount++;
                lnsCap[curArc] = cacheSize; 
                lnsCost[curArc] = 0;
                lnsToGArcId[curArc] = curEntry.innerArcId;
                //                LOG("curEA",i,curEntry.id,lnsG.id(curNode));
            }

        }

        for(auto it: lastSeen) {
            LOG("lastSeen",it.first.first,it.first.second,it.second);
        }
        assert(lastSeen.size()==0);

        LOG("extranode lastnode",lnsG.id(extraNode),lnsG.id(curNode),0);
        SmartDigraph::ArcIt aIt;
        /*        SmartDigraph::ArcIt aIt(lnsG);
        do {
            LOG("l arc",
                lnsToGNodeId[lnsG.source(aIt)],
                lnsToGNodeId[lnsG.target(aIt)],lnsCost[aIt]);
                } while (++aIt!=INVALID);*/



        LOG("createdLNS",nodeCount,arcCount,extraArcCount);
        double solval = solveMCF(lnsG, lnsCap, lnsCost, lnsSupplies, lnsFlow, 4, lnsPi);
        LOG("solLNS",solval,ejectNodes.size(),0);
        //        std::cerr << "ExLP" << 4 << " " << cacheSize << " hitc " << totalReqc-totalUniqC-solval << " reqc " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << std::endl;

        // recompute node potentials!
        
        // create residual graph
        SmartDigraph resG; // mcf graph
        SmartDigraph::NodeMap<int64_t> lnsToResNodeId(lnsG);
        SmartDigraph::NodeMap<int64_t> resToLnsNodeId(resG);
        SmartDigraph::ArcMap<double> resLength(resG);

        SmartDigraph::NodeIt nIt(lnsG);
        do {
            prevNode = curNode; // later needed to start bell solver
            curNode = resG.addNode();
            lnsToResNodeId[nIt] = resG.id(curNode);
            resToLnsNodeId[curNode] = lnsG.id(nIt);
        } while (++nIt!=INVALID);

        aIt = SmartDigraph::ArcIt(lnsG);
        do {
            bool addedArc = false;
            if(lnsFlow[aIt]<lnsCap[aIt]) {
                // add forward arc
                const int64_t resSourceId = lnsToResNodeId[lnsG.source(aIt)];
                const int64_t resTargetId = lnsToResNodeId[lnsG.target(aIt)];
                curArc = resG.addArc(resG.nodeFromId(resSourceId),
                                     resG.nodeFromId(resTargetId));
                resLength[curArc] = lnsCost[aIt];
                addedArc = true;
            }
            if(lnsFlow[aIt]>0) {
                // add reverse arc
                const int64_t resSourceId = lnsToResNodeId[lnsG.target(aIt)]; // reverse order
                const int64_t resTargetId = lnsToResNodeId[lnsG.source(aIt)];
                curArc = resG.addArc(resG.nodeFromId(resSourceId),
                                     resG.nodeFromId(resTargetId));
                resLength[curArc] = -lnsCost[aIt];
                addedArc = true;
            }
            assert(lnsCap[aIt]>0);
            assert(addedArc);
        } while (++aIt!=INVALID);

        BellSolveType bellsolver(resG, resLength);
        SmartDigraph::NodeMap<double> distResMap(resG);
        bellsolver.distMap(distResMap);

        /*        SmartDigraph::ArcIt raIt(resG);
        do {
            LOG("res arc",
                lnsToGNodeId[lnsG.nodeFromId(resToLnsNodeId[resG.source(raIt)])],
                lnsToGNodeId[lnsG.nodeFromId(resToLnsNodeId[resG.target(raIt)])],
                resLength[raIt]);
                } while (++raIt!=INVALID);*/

        
        LOG("bellsolver",lnsSupplies[lnsG.nodeFromId(resToLnsNodeId[prevNode])],0,0);
        bellsolver.run(prevNode);


        long double testVal;
        testVal = 0;
        // pi sum part
        SmartDigraph::NodeIt rnIt(resG);
        do {
            const int64_t lnsNodeId = resToLnsNodeId[rnIt];
            //            LOG("pi val supply",lnsNodeId,0,lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
            if(lnsSupplies[lnsG.nodeFromId(lnsNodeId)]!=0) {
                // ignore node that wouldn't count anyways
                // also conveniently ignores extraNode
                testVal += -distResMap[rnIt] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)];
                //                LOG("pi val calc",lnsNodeId,-distResMap[rnIt],lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
                if(!isfinite(testVal)) {
                    LOG("!!!nan",-distResMap[rnIt],lnsSupplies[lnsG.nodeFromId(lnsNodeId)],-distResMap[rnIt] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
                }
            }
        } while (++rnIt!=INVALID);
        LOG("dual val pi part",testVal,0,0);

        // alpha sum part
        SmartDigraph::ArcMap<double> lnsAlpha3(lnsG);
        aIt = SmartDigraph::ArcIt(lnsG);
        do {
            // ignore extraNode arc's alphas
            if(lnsCap[aIt] < maxFlowAmount) {
                const int64_t resSourceId = lnsToResNodeId[lnsG.source(aIt)];
                const int64_t resTargetId = lnsToResNodeId[lnsG.target(aIt)];
                const double piI = -distResMap[lnsG.nodeFromId(resSourceId)];
                const double piJ = -distResMap[lnsG.nodeFromId(resTargetId)];
                const double cij = lnsCost[aIt];
                const double curAlpha = piI - piJ - cij;
                if(curAlpha > 0)
                    lnsAlpha3[aIt] = curAlpha;
                else
                    lnsAlpha3[aIt] = 0;
                //                LOG("3rdalphas",lnsG.id(aIt),lnsAlpha3[aIt],lnsCap[aIt]);
            }
        } while (++aIt!=INVALID);

        long double testAlphaPart = 0;

        // iterate over lnsG arcs (as resG has extra arcs)
        aIt = SmartDigraph::ArcIt(lnsG);
        do {
            // ignore extra arc's alphas
            if(lnsCap[aIt] < maxFlowAmount) {
                testVal -= lnsAlpha3[aIt] * lnsCap[aIt];
                testAlphaPart -= lnsAlpha3[aIt] * lnsCap[aIt];
            }
        } while (++aIt!=INVALID);
        LOG("dual Val",testVal,testAlphaPart,totalReqc-totalUniqC-testVal);

        
        return 0;



    }


    return 0;
}

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
        LOG("g node",g.id(n),0,0);
    }
    std::cerr << "created graph with ";
    std::cerr << nodeCount << " nodes ";
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++arcCount;
        alpha[a] = 0;
        dualSol -= alpha[a] * cap[a];
        totalCapacity+=cap[a];
        const double curCost = cost[a];
        if(highestCost > curCost || highestCost==-1) {
            highestCost = curCost;
        }
        LOG("g arc",g.id(g.source(a)),g.id(g.target(a)),cost[a]);
    }
    std::cerr << arcCount << " arcs with dual cost " << dualSol << std::endl;

    std::cerr << std::endl << std::endl;

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
        std::vector<size_t> ejectSet;
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
        SmartDigraph::ArcMap<int64_t> lnsToGArcId(lnsG);

        // mark nodes for the ejection set
        traceIndex = kmin;
        arcId = -1;
        maxFlowAmount = 0;
        while(ejectSet.size()<maxEjectSize && traceIndex<trace.size()) {
            const trEntry & curEntry = trace[traceIndex];
            //            LOG("lnsG",traceIndex,curEntry.hasNext,trace.size());
            if(curEntry.hasNext) {
                ejectSet.push_back(traceIndex);
                LOG("ejectSet",traceIndex,0,0);
                arcId = curEntry.innerArcId;
                const int64_t ejNodeId = g.id(g.source(g.arcFromId(arcId)));
                ejectNodes.insert(ejNodeId);
                // calc max flow amount as sum of pos supplies
                maxFlowAmount += curEntry.size;
            }
            traceIndex++;
        }
        assert(maxFlowAmount>0);
        assert(arcId >= 0);
        LOG("ejS",kmin,maxEjectSize,ejectSet.size());

        // exception for last iteration
        if(traceIndex == trace.size()) {
            LOG("et",traceIndex,trace.size(),0);
            size_t i;
            for(i=traceIndex-1; !(trace[i].hasNext) && i>0;i--) {
            }
            arcId = trace[i].innerArcId;
            const int64_t ejNodeId1 = g.id(g.source(g.arcFromId(arcId)));
            const int64_t ejNodeId2 = g.id(g.target(g.arcFromId(arcId)));
            LOG("e node",ejNodeId1,ejNodeId2,0);
            ejectNodes.insert(ejNodeId2);
        }

        /*        for(auto it: ejectNodes) {
            LOG("lns g node",it,0,0);
            }*/

        LOG("ejS",kmin,maxEjectSize,ejectSet.size());

        // maybe iterate as in the parse script
        for(size_t i: ejectSet) {
            trEntry & curEntry = trace[i];
            auto curIdSize = std::make_pair(curEntry.id,curEntry.size);
            // reset extra arc flag
            extraArcFlag = false;
            extraArcCost = highestCost+1;

            // handle outer arc
            arcId = curEntry.outerArcId;

            // remember current node for outer arcid if in ejection set
            curArc = g.arcFromId(arcId);
            const int64_t iOutNodeId = g.id(g.target(curArc));
            LOG("iOutnodeid",iOutNodeId,0,arcId);
            if(ejectNodes.count(iOutNodeId)>0) {
                LOG("lastseen+iOutnodeid",iOutNodeId,0,arcId);
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


            // create all the original graph arcs

            // create outer arc
            const int64_t lNodeId = g.id(lNode);
            const uint64_t size = curEntry.size;
            if(lastSeen.count(curIdSize) > 0) {
                const SmartDigraph::Node lastReq = lnsG.nodeFromId(lastSeen[curIdSize]);
                curArc = lnsG.addArc(lastReq,curNode);
                arcCount++;
                lnsCap[curArc] = size;
                lnsCost[curArc] = 1/static_cast <double>(size);
                lnsSupplies[lastReq] += size;
                lnsSupplies[curNode] -= size;
                lnsToGArcId[curArc] = curEntry.outerArcId;
                lastSeen.erase(curIdSize);
            }

            // save mapping of curNode to g graph node
            lnsToGNodeId[curNode] = lNodeId;
            
            // inner arc
            prevNode = curNode;
            curNode = lnsG.addNode(); // next node
            nodeCount++;
            curArc = lnsG.addArc(prevNode,curNode);
            arcCount++;
            lnsCap[curArc] = cacheSize; 
            lnsCost[curArc] = 0;
            lnsToGArcId[curArc] = curEntry.innerArcId;

            // do we need to do second node??

        }

        LOG("are there no left lastSeens?",lastSeen.size(),0,0);
        
        // there will be outer arcs that remain to be added
        for(auto it: lastSeen) {
            const uint64_t size = it.first.second;
            LOG("adding lastSeen remainder",it.first.first,size,it.second);
            const SmartDigraph::Node lastReq = lnsG.nodeFromId(it.second);
            curArc = lnsG.addArc(lastReq,curNode);
            arcCount++;
            lnsCap[curArc] = size;
            lnsCost[curArc] = 1/static_cast <double>(size);
            lnsSupplies[lastReq] += size;
            lnsSupplies[curNode] -= size;
            // BUG FIXME
            //lnsToGArcId[curArc] = curEntry.outerArcId;
        }


        LOG("extranode lastnode",lnsG.id(extraNode),lnsG.id(curNode),0);
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            LOG("l arcs",
                lnsToGNodeId[lnsG.source(it)],
                lnsToGNodeId[lnsG.target(it)],lnsCost[it]);
        }

        LOG("createdLNS",nodeCount,arcCount,extraArcCount);
        double solval = solveMCF(lnsG, lnsCap, lnsCost, lnsSupplies, lnsFlow, 4, lnsPi);
        LOG("solLNS",solval,ejectNodes.size(),0);
        //        std::cerr << "ExLP" << 4 << " " << cacheSize << " hitc " << totalReqc-totalUniqC-solval << " reqc " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << std::endl;

        long double testVal;

        // create residual graph
        SmartDigraph resG; // mcf graph
        SmartDigraph::NodeMap<int64_t> lnsResMap(lnsG);
        SmartDigraph::NodeMap<int64_t> resLnsMap(resG);
        SmartDigraph::ArcMap<double> resLength(resG);
        for (SmartDigraph::NodeIt it(lnsG); it!=INVALID; ++it) {
            curNode = resG.addNode();
            lnsResMap[it] = resG.id(curNode);
            resLnsMap[curNode] = lnsG.id(it);
        }

        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            bool addedArc = false;
            if(lnsFlow[it]<lnsCap[it]) {
                // add forward arc
                const int64_t resSourceId = lnsResMap[lnsG.source(it)];
                const int64_t resTargetId = lnsResMap[lnsG.target(it)];
                curArc = resG.addArc(resG.nodeFromId(resSourceId),
                                     resG.nodeFromId(resTargetId));
                resLength[curArc] = lnsCost[it];
                addedArc = true;
            }
            if(lnsFlow[it]>0) {
                // add reverse arc
                const int64_t resSourceId = lnsResMap[lnsG.target(it)]; // reverse order
                const int64_t resTargetId = lnsResMap[lnsG.source(it)];
                curArc = resG.addArc(resG.nodeFromId(resSourceId),
                                     resG.nodeFromId(resTargetId));
                resLength[curArc] = -lnsCost[it];
                addedArc = true;
            }
            assert(lnsCap[it]>0);
            assert(addedArc);
        }

        /*for (SmartDigraph::ArcIt it(resG); it!=INVALID; ++it) {
            //            LOG("resG",resG.id(resG.source(it)),resG.id(resG.target(it)),lnsCost[it]);//resG.id(it),
            LOG("resG",resLnsMap[resG.source(it)],resLnsMap[resG.target(it)],resLength[it]);
            }*/
        
        BellSolveType bellsolver(resG, resLength);
        SmartDigraph::NodeMap<double> distResMap(resG);
        bellsolver.distMap(distResMap);
        
        bellsolver.run(curNode);

        /*        for (SmartDigraph::NodeIt it(resG); it!=INVALID; ++it) {
            LOG("belldists",resG.id(it),-distResMap[it],0);
        }
        */

        testVal = 0;
        // pi sum part
        for (SmartDigraph::NodeIt it(resG); it!=INVALID; ++it) {
            const int64_t lnsNodeId = resLnsMap[it];
            testVal += -distResMap[it] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)];
            //LOG("3rdval",lnsNodeId,-distResMap[it],lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
            if(!isfinite(testVal)) {
                LOG("!!!nan",-distResMap[it],lnsSupplies[lnsG.nodeFromId(lnsNodeId)],-distResMap[it] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
            }

        }
        // alpha sum part
        SmartDigraph::ArcMap<double> lnsAlpha3(lnsG);
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            const int64_t resSourceId = lnsResMap[lnsG.source(it)];
            const int64_t resTargetId = lnsResMap[lnsG.target(it)];
            const double piI = -distResMap[lnsG.nodeFromId(resSourceId)];
            const double piJ = -distResMap[lnsG.nodeFromId(resTargetId)];
            const double cij = lnsCost[it];
            const double curAlpha = piI - piJ - cij;
            if(curAlpha > 0)
                lnsAlpha3[it] = curAlpha;
            else
                lnsAlpha3[it] = 0;
            //LOG("3rdalphas",lnsG.id(it),lnsAlpha3[it],lnsCap[it]);
        }

        testVal = 0;
        long double testAlphaPart = 0, testAlphaPartAll = 0;
        // iterate over resG nodes (as they are the same as lnsG nodes)
        for (SmartDigraph::NodeIt it(resG); it!=INVALID; ++it) {
            const int64_t lnsNodeId = resLnsMap[it];
            testVal += -distResMap[it] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)];
        }

        // iterate over lnsG arcs (as resG has extra arcs)
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            // ignore extra arc's alphas
            if(lnsCap[it] < maxFlowAmount) {
                testVal -= lnsAlpha3[it] * lnsCap[it];
                testAlphaPart -= lnsAlpha3[it] * lnsCap[it];
            }
            testAlphaPartAll -= lnsAlpha3[it] * lnsCap[it];
        }
        LOG("dual Val",testVal,testAlphaPart,testAlphaPartAll);

        
        return 0;



    }


    return 0;
}

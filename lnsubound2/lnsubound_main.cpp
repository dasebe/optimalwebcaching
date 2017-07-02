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

// uncomment to enable cache debugging:
#define CDEBUG 1
// util for debug
#ifdef CDEBUG
inline void logMessage(std::string m, double x, double y, double z) {
    std::cerr << m << "\t" << x << "\t" << y  << "\t" << z << "\n";
}
#define LOG(m,x,y,z) logMessage(m,x,y,z)
#else
#define LOG(m,x,y,z)
#endif

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
    SmartDigraph::Node lastNode = createMCF(g, trace, cacheSize, cap, cost, supplies);
    
    // LNS constants
    int64_t totalCapacity = 0;
    int64_t maxFlowAmount = 0;
    double highestCost = -1;

    // dual-feasible solution: all alpha and pi are 0
    uint64_t nodes = 0, vertices = 0;
    double dualSol = 0;
    SmartDigraph::NodeMap<double> pi(g);
    SmartDigraph::ArcMap<double> alpha(g);

    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodes;
        pi[n] = 0;
        dualSol += pi[n] * supplies[n];
    }
    std::cerr << "created graph with ";
    std::cerr << nodes << " nodes ";
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++vertices;
        alpha[a] = 0;
        dualSol -= alpha[a] * cap[a];
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
    
    // LNS iteration step state
    std::vector<size_t> ejectSet;
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
    size_t traceIndex;

    // todo fixme: kmin
    for(size_t kmin=4; kmin<trace.size(); kmin=traceIndex - maxEjectSize/2) {
        // LNS graph structure
        SmartDigraph lnsG; // mcf graph
        extraNode = lnsG.addNode();
        curNode = lnsG.addNode();
        SmartDigraph::ArcMap<int64_t> lnsCap(lnsG); // mcf capacities
        SmartDigraph::ArcMap<double> lnsCost(lnsG); // mcf costs
        SmartDigraph::NodeMap<int64_t> lnsSupplies(lnsG); // mcf demands/supplies
        SmartDigraph::ArcMap<int64_t> lnsFlow(lnsG); // in units of flow
        SmartDigraph::NodeMap<double> lnsPi(lnsG);
        SmartDigraph::ArcMap<double> lnsAlpha(lnsG);

        // mark nodes for the ejection set
        traceIndex = kmin;
        arcId = -1;
        maxFlowAmount = 0;
        while(ejectSet.size()<maxEjectSize && traceIndex<trace.size()) {
            const trEntry & curEntry = trace[traceIndex];
            if(curEntry.hasNext) {
                arcId = curEntry.innerArcId;
                const int64_t ejNodeId = g.id(g.source(g.arcFromId(arcId)));
                ejectSet.push_back(traceIndex);
                ejectNodes.insert(ejNodeId);
                // calc max flow amount as sum of pos supplies
                maxFlowAmount += curEntry.size;
            }
            traceIndex++;
        }
        LOG("ejS",kmin,maxEjectSize,ejectSet.size());
        assert(maxFlowAmount>0);
        assert(arcId >= 0);

        for(size_t i: ejectSet) {
            trEntry & curEntry = trace[i];
            // reset extra arc flag
            extraArcFlag = false;
            extraArcCost = highestCost+1;

            // handle outer arc
            arcId = curEntry.outerArcId;

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
                lnsCost[curArc] = extraArcCost;
                lnsCap[curArc] = maxFlowAmount;
            }


            // create all the original graph arcs

            // create outer arc
            const int64_t lNodeId = g.id(lNode);
            const uint64_t size = curEntry.size;
            if(lastSeen.count(lNodeId) > 0) {
                const SmartDigraph::Node lastReq = lnsG.nodeFromId(lastSeen[lNodeId]);
                curArc = lnsG.addArc(lastReq,curNode);
                lnsCap[curArc] = size;
                lnsCost[curArc] = 1/static_cast <double>(size);
                lnsSupplies[lastReq] += size;
                lnsSupplies[curNode] -= size;
            }

            // inner arc
            prevNode = curNode;
            curNode = lnsG.addNode(); // next node
            curArc = lnsG.addArc(prevNode,curNode);
            lnsCap[curArc] = cacheSize; 
            lnsCost[curArc] = 0;

            // do we need to do second node??

        }

        /*LOG("extranode lastnode",lnsG.id(extraNode),lnsG.id(curNode),0);
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            LOG("lnsG",lnsG.id(lnsG.source(it)),lnsG.id(lnsG.target(it)),lnsCost[it]);//lnsG.id(it),
            }*/


        double solval = solveMCF(lnsG, lnsCap, lnsCost, lnsSupplies, lnsFlow, 4, lnsPi);
        LOG("solLNS",solval,ejectNodes.size(),0);
        //        std::cerr << "ExLP" << 4 << " " << cacheSize << " hitc " << totalReqc-totalUniqC-solval << " reqc " << totalReqc << " OHR " << 1.0-(static_cast<double>(solval)+totalUniqC)/totalReqc << std::endl;

        /*        LOG("-----------nodes---------",0,0,0);
        for (SmartDigraph::NodeIt i(lnsG); i!=INVALID; ++i) {
            LOG("",lnsG.id(i),lnsSupplies[i],lnsPi[i]);
        }*/

        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            const double piI = lnsPi[lnsG.source(it)];
            const double piJ = lnsPi[lnsG.target(it)];
            const double cij = lnsCost[it];
            const double curAlpha = piI - piJ - cij;
            if(curAlpha > 0)
                lnsAlpha[it] = curAlpha;
            else
                lnsAlpha[it] = 0;
            //LOG("alphaC\t"+std::to_string(lnsAlpha[it]),lnsPi[lnsG.source(it)],lnsPi[lnsG.target(it)],lnsCost[it]);
        }

        long double testVal;


        testVal = 0;
        for (SmartDigraph::NodeIt it(lnsG); it!=INVALID; ++it) {
            testVal += lnsPi[it] * lnsSupplies[it];
            if(!isfinite(testVal)) {
                LOG("!!!nan",lnsPi[it],lnsSupplies[it],lnsPi[it] * lnsSupplies[it]);
            }
        }
        LOG("not quite dual Val",testVal,0,0);

        testVal = 0;
        for (SmartDigraph::NodeIt it(lnsG); it!=INVALID; ++it) {
            testVal += lnsPi[it] * lnsSupplies[it];
            //            LOG("pis",lnsG.id(it),lnsPi[it],0);
        }
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            testVal -= lnsAlpha[it] * lnsCap[it];
            //LOG("alphas",lnsG.id(it),lnsAlpha[it],lnsCap[it]);
        }
        LOG("dual Val",testVal,0,0);

        testVal = 0;
        long double flowInCache = 0, flowOutCache = 0, extraFlow = 0;
        uint64_t arcIn = 0, arcOut = 0, arcExtraTo = 0, arcExtraFrom = 0;
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            testVal += lnsFlow[it] * lnsCost[it];
            if(lnsCap[it] == cacheSize) {
                flowInCache += lnsFlow[it]* lnsCost[it];
                arcIn++;
            } else if (lnsG.source(it) == extraNode) {
                extraFlow += lnsFlow[it]* lnsCost[it];
                arcExtraFrom++;
            } else if (lnsG.target(it) == extraNode) {
                extraFlow += lnsFlow[it]* lnsCost[it];
                arcExtraTo++;
            } else {
                flowOutCache += lnsFlow[it]* lnsCost[it];
                arcOut++;
            }
            //LOG("primal calc: flow, cost, cap ",lnsFlow[it],lnsCost[it],lnsCap[it]);
        }
        LOG("primal Val",testVal,0,0);
        LOG("primal cost: in,out,extra ",flowInCache,flowOutCache,extraFlow);
        LOG("arc counts: in,out,extrato,extrafrom\t"+std::to_string(arcIn),arcOut,arcExtraTo,arcExtraFrom);



        //// test

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
        for (SmartDigraph::NodeIt it(resG); it!=INVALID; ++it) {
            const int64_t lnsNodeId = resLnsMap[it];
            testVal += -distResMap[it] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)];
            //LOG("3rdval",lnsNodeId,-distResMap[it],lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
            if(!isfinite(testVal)) {
                LOG("!!!nan",-distResMap[it],lnsSupplies[lnsG.nodeFromId(lnsNodeId)],-distResMap[it] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)]);
            }

        }
        LOG("3rd not quite dual Val",testVal,0,0);
        

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
        for (SmartDigraph::NodeIt it(resG); it!=INVALID; ++it) {
            const int64_t lnsNodeId = resLnsMap[it];
            testVal += -distResMap[it] * lnsSupplies[lnsG.nodeFromId(lnsNodeId)];
        }

        // no need to check if extra node contributes to potential sum
        //        const int64_t eNodeId = lnsResMap[extraNode];
        //        LOG("check check",-distResMap[resG.nodeFromId(eNodeId)],lnsSupplies[extraNode],-distResMap[resG.nodeFromId(eNodeId)] * lnsSupplies[extraNode]);
        
        for (SmartDigraph::ArcIt it(lnsG); it!=INVALID; ++it) {
            //            LOG("2ndalphasl",lnsG.id(it),lnsAlpha3[it],lnsCap[it]);
            if(lnsCap[it] < maxFlowAmount) {
                //                LOG("2ndalphasc",lnsG.id(it),lnsAlpha3[it],lnsCap[it]);
                testVal -= lnsAlpha3[it] * lnsCap[it];
                testAlphaPart -= lnsAlpha3[it] * lnsCap[it];
            }
            testAlphaPartAll -= lnsAlpha3[it] * lnsCap[it];
        }
        LOG("3rd dual Val",testVal,testAlphaPart,testAlphaPartAll);

        
        //// test


        return 0;



    }


    return 0;
}

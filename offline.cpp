#include <iostream>
#include <cassert>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <lemon/smart_graph.h>
#include <lemon/lgf_writer.h>
#include <lemon/capacity_scaling.h>

using namespace lemon;

void createMCF(SmartDigraph & g, std::vector<std::tuple<uint64_t,uint64_t,bool> > & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int> & supplies) {

    // init
    std::unordered_map<uint64_t, uint64_t> lastSeen;
    uint64_t reqc=0;
    SmartDigraph::Arc curArc;
    SmartDigraph::Node curNode = g.addNode(); // initial node
    SmartDigraph::Node prevNode;
    std::unordered_map<uint64_t, uint64_t> lastSize;

    // iterate over trace
    for(auto it: trace) {
        // only consider requests that reoccur
        // should we delete if does not reoccur?
        if(std::get<2>(it)) {
            const uint64_t id=std::get<0>(it);
            const uint64_t size=std::get<1>(it);
            if(lastSeen.count(id)>0) {
                if(size==lastSize[id]) {
                    // create request arc ("outer")
                    const SmartDigraph::Node lastReq = g.nodeFromId(lastSeen[id]);
                    curArc = g.addArc(lastReq,curNode);
                    cap[curArc] = size;
                    cost[curArc] = 1/static_cast <double>(size);
                    supplies[lastReq] += size;
                    supplies[curNode] -= size;
                } else {
                    std::cerr << "fuckup" << std::endl;
                    assert(false);
                }
            }
            // create capacity arc ("inner")
            prevNode = curNode;
            curNode = g.addNode(); // next node
            curArc = g.addArc(prevNode,curNode);
            cap[curArc] = cacheSize; 
            cost[curArc] = 0;
            lastSeen[id]=reqc;
            lastSize[id]=size;
            reqc++;
        }
    }
    
    // there will m intervals without arcs
    // they all intersect and thus point to the same final node
    for(auto it: lastSeen) {
        const uint64_t id= it.second;
        const uint64_t size=lastSize[id];
        const SmartDigraph::Node lastReq = g.nodeFromId(id);
        curArc = g.addArc(lastReq,curNode);
        cap[curArc] = size;
        cost[curArc] = 1/static_cast <double>(size);
        supplies[lastReq] += size;
        supplies[curNode] -= size;
    }
}    

int main(int argc, char* argv[]) {

    std::string path(argv[1]);
    uint64_t cacheSize(atoll(argv[2]));
    std::ifstream traceFile(path);

    uint64_t t, id, size, reqc=0, uniqc=0;
    std::vector<std::tuple<uint64_t,uint64_t,bool> > trace;
    std::unordered_map<uint64_t, uint64_t> lastSeen;

    while(traceFile >> t >> id >> size) {
        if(lastSeen.count(id)>0) {
            std::get<2>(trace[lastSeen[id]]) = true;
        } else {
            uniqc++;
        }
        trace.emplace_back(id,size,false);
        lastSeen[id]=reqc++;
    }
    traceFile.close();
    lastSeen.clear();

    std::cout << "scanned trace n=" << reqc << " m= " << uniqc << std::endl;

    SmartDigraph g; // mcf graph
    g.reserveNode(reqc-uniqc+1);
    g.reserveArc(2*(reqc-uniqc));
    SmartDigraph::ArcMap<int> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int> supplies(g); // mcf demands/supplies
    createMCF(g, trace, cacheSize, cap, cost, supplies);
    trace.clear();
    
    std::cout << "created graph with ";

    int nodes=0, vertices=0;
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) ++nodes;
    std::cout << nodes << " nodes ";
    
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) ++vertices;
    std::cout << vertices << " arcs ";

    digraphWriter(g, std::cout).
        arcMap("capacity", cap).       // write cap into 'capacity'
        arcMap("cost", cost).          // write 'cost' for for arcs
        run();
    
    return 0;
}

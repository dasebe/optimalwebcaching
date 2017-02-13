#include <string>
#include <list>
#include <lemon/lgf_writer.h>
#include "catch.hpp"
#include "../offline.h"

TEST_CASE( "MCF create function", "[offline]" ) {

    std::vector<std::tuple<uint64_t,uint64_t,bool> > trace;
    std::string path("test.tr");
    parseTraceFile(trace, path);
    REQUIRE(trace.size()==8);

    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int> supplies(g); // mcf demands/supplies
    uint64_t cacheSize = 10;
    createMCF(g, trace, cacheSize, cap, cost, supplies);

    digraphWriter(g, std::cout).
        nodeMap("supplies", supplies).
        arcMap("capacity", cap).       // write cap into 'capacity'
        arcMap("cost", cost).          // write 'cost' for for arcs
        run();

    int nodes=0, vertices=0, supplysum=0;
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodes;
        supplysum+=supplies[n];
    }
    REQUIRE(nodes==6);
    REQUIRE(supplysum==0);
    
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++vertices;
    }
    REQUIRE(vertices==10);

    // check details
    // Checking supply/demand values
    REQUIRE(supplies[g.nodeFromId(0)]==2);
    REQUIRE(supplies[g.nodeFromId(1)]==3);
    REQUIRE(supplies[g.nodeFromId(2)]==0);
    REQUIRE(supplies[g.nodeFromId(3)]==4);
    REQUIRE(supplies[g.nodeFromId(4)]==0);
    REQUIRE(supplies[g.nodeFromId(5)]==-1*(2+3+4));

    // Checking costs and capacities
    std::list<uint64_t> cl {0, 1, 3, 4, 6};
    for(auto it: cl) {
        REQUIRE(cap[g.arcFromId(it)]==10);
        REQUIRE(cost[g.arcFromId(it)]==0);
    }
    std::list<uint64_t> id1 {2, 5, 7};
    for(auto it: id1) {
        REQUIRE(cap[g.arcFromId(it)]==2);
        REQUIRE(cost[g.arcFromId(it)]==1/2.0);
    }
    REQUIRE(cap[g.arcFromId(8)]==3);
    REQUIRE(cost[g.arcFromId(8)]==1/3.0);
    REQUIRE(cap[g.arcFromId(9)]==4);
    REQUIRE(cost[g.arcFromId(9)]==1/4.0);
}

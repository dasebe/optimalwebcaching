#include <string>
#include <list>
#include <lemon/lgf_writer.h>
#include "catch.hpp"
#include "../lib/parse_trace.h"

TEST_CASE( "test trace 1: basic MCF graph","[trace1]") {

    std::vector<traceEntry> trace;
    std::string path("test1.tr");
    uint64_t uniqc = parseTraceFile(trace, path);
    REQUIRE(uniqc==2);
    REQUIRE(trace.size()==4);

    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    uint64_t cacheSize = 2;
    createMCF(g, trace, cacheSize, cap, cost, supplies);

    std::cout << "---------test1 graph--------" << std::endl;
    digraphWriter(g, std::cout).
        nodeMap("supplies", supplies).
        arcMap("capacity", cap).       // write cap into 'capacity'
        arcMap("cost", cost).          // write 'cost' for for arcs
        run();

    uint64_t nodes=0, vertices=0, supplysum=0;
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodes;
        supplysum+=supplies[n];
    }
    REQUIRE(nodes==3);
    REQUIRE(supplysum==0);
    
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++vertices;
    }
    REQUIRE(vertices==4);
    
    std::list<uint64_t> cl {0, 2};
    for(auto it: cl) {
        REQUIRE(cap[g.arcFromId(it)]==2);
        REQUIRE(cost[g.arcFromId(it)]==0);
    }

    REQUIRE(cap[g.arcFromId(1)]==2);
    REQUIRE(cost[g.arcFromId(1)]==1/2.0);
    REQUIRE(cap[g.arcFromId(3)]==3);
    REQUIRE(cost[g.arcFromId(3)]==1/3.0);
}


TEST_CASE( "test trace 2: larger MCF graph","[trace2]") {

    std::vector<traceEntry> trace;
    std::string path("test2.tr");
    uint64_t uniqc = parseTraceFile(trace, path);
    REQUIRE(uniqc==3);
    REQUIRE(trace.size()==8);

    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    uint64_t cacheSize = 10;
    createMCF(g, trace, cacheSize, cap, cost, supplies);

    std::cout << "---------test2 graph--------" << std::endl;
    digraphWriter(g, std::cout).
        nodeMap("supplies", supplies).
        arcMap("capacity", cap).       // write cap into 'capacity'
        arcMap("cost", cost).          // write 'cost' for for arcs
        run();

    uint64_t nodes=0, vertices=0, supplysum=0;
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
    std::list<uint64_t> id1 {2, 5, 8};
    for(auto it: id1) {
        REQUIRE(cap[g.arcFromId(it)]==2);
        REQUIRE(cost[g.arcFromId(it)]==1/2.0);
    }
    REQUIRE(cap[g.arcFromId(7)]==3);
    REQUIRE(cost[g.arcFromId(7)]==1/3.0);
    REQUIRE(cap[g.arcFromId(9)]==4);
    REQUIRE(cost[g.arcFromId(9)]==1/4.0);
}





TEST_CASE( "test trace 3: MCF graph with id/size inconsistency","[trace3]") {

    std::vector<traceEntry> trace;
    std::string path("test3.tr");
    uint64_t uniqc = parseTraceFile(trace, path);
    REQUIRE(uniqc==13); //12 ids and one size inconsistency
    REQUIRE(trace.size()==15);

    SmartDigraph g; // mcf graph
    SmartDigraph::ArcMap<int64_t> cap(g); // mcf capacities
    SmartDigraph::ArcMap<double> cost(g); // mcf costs
    SmartDigraph::NodeMap<int64_t> supplies(g); // mcf demands/supplies

    uint64_t cacheSize = 2;
    createMCF(g, trace, cacheSize, cap, cost, supplies);

    std::cout << "---------test3 graph--------" << std::endl;
    digraphWriter(g, std::cout).
        nodeMap("supplies", supplies).
        arcMap("capacity", cap).       // write cap into 'capacity'
        arcMap("cost", cost).          // write 'cost' for for arcs
        run();

    uint64_t nodes=0, vertices=0, supplysum=0;
    for (SmartDigraph::NodeIt n(g); n!=INVALID; ++n) {
        ++nodes;
        supplysum+=supplies[n];
    }
    REQUIRE(nodes==3);
    REQUIRE(supplysum==0);
    
    for (SmartDigraph::ArcIt a(g); a != INVALID; ++a) {
        ++vertices;
    }
    REQUIRE(vertices==4);
    
    std::list<uint64_t> cl {0, 2};
    for(auto it: cl) {
        REQUIRE(cap[g.arcFromId(it)]==2);
        REQUIRE(cost[g.arcFromId(it)]==0);
    }

    REQUIRE(cap[g.arcFromId(1)]==4294967297);
    REQUIRE(cost[g.arcFromId(1)]==1/4294967297.0);
    REQUIRE(cap[g.arcFromId(3)]==1);
    REQUIRE(cost[g.arcFromId(3)]==1/1.0);
}

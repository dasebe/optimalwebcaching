#include <cassert>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <math.h>
#include "lib/parse_trace.h"

int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << argv[0] << " traceFile" << std::endl;
        return 1;
    }

    std::string path(argv[1]);

    // parse trace file
    std::vector<trEntry> trace;
    parseTraceFile(trace, path);
    //    uint64_t totalReqc = trace.size();

    std::map<double, uint64_t> uniqsizedist;

    // get nextSeen indices
    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t> lastSeen;
    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t> reqcounter;
    for (size_t i = trace.size(); i--> 0 ;) {
        //for(int64_t i=trace.size()-1; i>=0; i--) {
        trEntry & cur = trace[i];
        if(lastSeen.count(std::make_pair(cur.id,cur.size)) > 0) {
            cur.hasNext = true;
            cur.nextSeen = lastSeen[std::make_pair(cur.id, cur.size)];
        } else {
            uniqsizedist[std::round(std::log10(cur.size)*10.0)/10.0]++;
        }
        lastSeen[std::make_pair(cur.id, cur.size)] = i;
        reqcounter[std::make_pair(cur.id, cur.size)]++;
    }

    // overall reuse distance
    std::map<double, uint64_t> rddist;
    for(size_t i=0; i<trace.size(); i++) {
        auto ns = trace[i].nextSeen;
        if(ns==0) {
            rddist[-1.0]++;
        } else {
            rddist[std::round(std::log10(ns-i)*100.0)/100.0]++;
        }
    }

    for(auto it:rddist) {
        std::cout << "rd 0 " << it.first << " " << it.second << "\n";
    }

    // overall popularity
    std::map<double, uint64_t> popdist;
    std::map<double, uint64_t> sizedist;
    for(auto it:reqcounter) {
        popdist[std::round(std::log10(it.second)*100.0)/100.0]++;
        sizedist[std::round(std::log10(std::get<1>(it.first))*100.0)/100.0]+=it.second;
    }

    for(auto it:popdist) {
        std::cout << "pop 0 " << it.first << " " << it.second << "\n";
    }
        
    for(auto it:sizedist) {
        std::cout << "size 0 " << it.first << " " << it.second << "\n";
    }

    for(auto it:uniqsizedist) {
        std::cout << "uniqsize 0 " << it.first << " " << it.second << "\n";
    }


    std::vector<uint64_t> zipfdist;
    for(auto it:reqcounter) {
        zipfdist.push_back(it.second);
    }
    std::sort(zipfdist.begin(),zipfdist.end());
    size_t pos=1, printpos=1;
    for(auto it=zipfdist.crbegin();it!=zipfdist.crend();it++){
        if(pos==printpos) {
            std::cout << "zipf 0 " << pos << " " << *it << "\n";
            printpos*=2;
        }
        pos++;
    }
    return 0;
}



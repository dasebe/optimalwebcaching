#include <iostream>
#include <random>
#include <fstream>
#include <cassert>
#include <map>
#include <vector>
#include <unordered_map>
#include <set>
#include <tuple>
#include "parse_trace.h"

void parseTraceFile(std::string & path, std::string & pathout, uint64_t samples) {
    std::ifstream traceFile(path);
   
    uint64_t time, id, size, reqc=0;

    std::map<std::pair<uint64_t, uint64_t>, uint64_t> counts;

    while(traceFile >> time >> id >> size) {
        counts[std::make_pair(id,size)]++;
	reqc++;
    }
    traceFile.clear();
    traceFile.seekg(0, std::ios::beg);

    std::default_random_engine generator;

    std::vector<std::pair<uint64_t, uint64_t> > ids;
    for(auto it: counts) {
      ids.push_back(it.first);
    }

    std::uniform_int_distribution<size_t> distribution(0,ids.size()-1);
    if(reqc < samples/2) {
      std::cerr << "error: too many samples\n";
      assert(reqc> samples/2);
    }

    std::set<std::pair<uint64_t, uint64_t>> sampled;
    uint64_t sampleCount = 0;

    while(true) {
      const auto idx = distribution(generator);
      const auto idSize = ids[idx];
      if(sampleCount + counts[idSize] > samples) {
	break;
      }
      sampleCount += counts[idSize];
      sampled.insert(idSize);
    }

    std::ofstream outFile(pathout);

    while(traceFile >> time >> id >> size) {
        auto idSize = std::make_pair(id,size);
        // skip non-sampled objects
        if(sampled.count(idSize)==0) {
            continue;
        }
	outFile << time << " " << id << " " << size << "\n";

    }
}

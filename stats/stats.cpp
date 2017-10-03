#include <cassert>
#include <vector>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include "lib/parse_trace.h"


int main(int argc, char* argv[]) {

    if (argc != 2) {
        std::cerr << argv[0] << " traceFile" << std::endl;
        return 1;
    }

    std::cerr << "starting..\n";

    std::string path(argv[1]);
    std::vector<trEntry> trace;
    std::ifstream traceFile(path);
    uint64_t time, id, size, reqc=0;

    std::cerr << "parsing trace\n";

    std::unordered_map<std::pair<uint64_t, uint64_t>, uint64_t> lastSeen;
    while(traceFile >> time >> id >> size) {
        const auto idSize = std::make_pair(id,size);
        if(lastSeen.count(idSize)>0) {
            trace[lastSeen[idSize]].nextSeen = reqc;
        }
        trace.emplace_back(size);
        lastSeen[idSize]=reqc++;
    }

    std::cerr << "parsing intervals\n";

    std::unordered_map<size_t, size_t> interval;
    for(size_t i=0; i<trace.size(); i++) {
      trEntry & cur = trace[i];
      
      // check if I have any parents
      if(cur.nextSeen > 0) {
	for(auto it: interval) {
	  if( (it.first > cur.nextSeen) && (trace[it.second].size >= cur.size) ) {
	    // this is my parent
	    cur.parents.insert(it.second);
	  }
	}
      }

      if(interval.count(i)>0) {
	// end this interval
	interval.erase(i);
      }

      if(cur.nextSeen > 0) {
	// start new interval
	interval[cur.nextSeen] = i;
      }
    }

    std::cerr << "outputting\n";

    for(size_t i=0; i<trace.size(); i++) {
      for(size_t parent: trace[i].parents) {
	std::cout << i << " " << trace[i].nextSeen << " " << trace[i].size << " - " << parent << " " << trace[parent].nextSeen << " " << trace[parent].size <<"\n";
      }
    }

    return 0;
}

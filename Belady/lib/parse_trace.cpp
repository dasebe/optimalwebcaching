//#include <iostream>
#include <fstream>
#include <map>
#include <cassert>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include "parse_trace.h"

void parseTraceFile(std::vector<trEntry> & trace, std::string & path) {
    std::ifstream traceFile(path);
    uint64_t time, id, size;
    while(traceFile >> time >> id >> size) {
        trace.emplace_back(id,size);
    }
}

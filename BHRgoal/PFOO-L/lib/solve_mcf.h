#include "parse_trace.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

void cacheAlg(std::vector<trEntry> & trace);
void printRes(std::vector<trEntry> & trace, uint64_t & byteSum, uint64_t cacheSize);

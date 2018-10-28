#include "parse_trace.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

void cacheAlg(std::vector<trEntry> & trace);
void printRes(std::vector<trEntry> & trace, std::string algName, uint64_t cacheSize, std::ofstream * r);

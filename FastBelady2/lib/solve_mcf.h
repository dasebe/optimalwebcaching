#include "parse_trace.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>

void cacheAlg(std::vector<trEntry*> & trace, uint64_t cacheSize);
void printRes(std::vector<trEntry*> & trace, std::string algName);

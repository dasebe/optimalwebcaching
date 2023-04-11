#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <fstream>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

int main(int argc, char* argv[]) {

    if (argc != 5) {
        std::cerr << argv[0] << " traceFile cacheSizeMax resultPath warmup" << std::endl;
        return 1;
    }

    std::cerr << "start up\n";
    std::string path(argv[1]);
    uint64_t cacheSizeMax(std::stoull(argv[2]));
    std::string resultPath(argv[3]);
    uint64_t warmup(std::stoull(argv[4]));
    std::ofstream * resultFile = new std::ofstream(resultPath);

    // parse trace file
    std::vector<trEntry> trace;
    parseTraceFile(trace, path, warmup);
    std::cerr << "parsed up\n";

    cacheAlg(trace);

    std::cerr << "sorted up\n";
    
    printRes(trace, "fluid2", cacheSizeMax, resultFile);

    return 0;
}



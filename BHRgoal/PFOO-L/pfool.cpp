#include <cassert>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <fstream>
#include "lib/parse_trace.h"
#include "lib/solve_mcf.h"

int main(int argc, char* argv[]) {

    if (argc != 3) {
        std::cerr << argv[0] << " traceFile cacheSize" << std::endl;
        return 1;
    }

    std::cerr << "start up\n";
    std::string path(argv[1]);
    uint64_t cacheSize(std::stoull(argv[2]));

    // parse trace file
    std::vector<trEntry> trace;
    uint64_t byteSum;
    parseTraceFile(trace, path,byteSum);
    std::cerr << "parsed " << byteSum << " bytes\n";

    cacheAlg(trace);
    std::cerr << "processed\n";
    
    printRes(trace, byteSum, cacheSize);

    return 0;
}



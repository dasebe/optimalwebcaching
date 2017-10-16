#include <fstream>
#include <iostream>
#include <cassert>
#include "lib/parse_trace.h"

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::cerr << argv[0] << " traceFile outFile samples" << std::endl;
        return 1;
    }


    std::string path(argv[1]);
    std::cerr << "starting " << path << "\n";
    std::string outPath(argv[2]);
    uint64_t samples(atoll(argv[3]));

    parseTraceFile(path, outPath, samples);

    std::cerr << "done " << path << "\n";

    return 0;
}

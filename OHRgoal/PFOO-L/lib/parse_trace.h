#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <misc/hash_combine.h>

// uncomment to enable debugging:
//#define DEBUG 1
#ifdef DEBUG
#define LOG(m,x,y,z) log_message(m,x,y,z,"\n")
#else
#define LOG(m,x,y,z)
#endif
inline void log_message(std::string m, double x, double y, double z, std::string e) {
    std::cerr << m << "," << x << "," << y  << "," << z << e;
}

// trace entry
struct trEntry {
    uint64_t size;
    uint64_t volume;
    bool hasNext;

    trEntry(uint64_t nsize)
        : size(nsize),
          volume(std::numeric_limits<uint64_t>::max()),
          hasNext(false)
    {
    };

    bool operator <(const trEntry &b) //const trEntry &a
    {
        return volume < b.volume;
    };
};

void parseTraceFile(std::vector<trEntry> & trace, std::string & path, uint64_t warmup);

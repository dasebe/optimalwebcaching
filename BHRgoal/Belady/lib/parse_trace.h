#pragma once
#include <misc/hash_combine.h>
#include <vector>
#include <string>
#include <iostream>

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
    const uint64_t id;
    const uint64_t size;
    size_t nextSeen;
    bool hasNext;
    bool hit;

    trEntry(uint64_t nid, uint64_t nsize)
        : id(nid),
          size(nsize),
          nextSeen(0),
          hasNext(false),
          hit(false)
    {
    };
};

void parseTraceFile(std::vector<trEntry> & trace, std::string & path);

#pragma once
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

// from boost hash combine: hashing of std::pairs for unordered_maps
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}

void parseTraceFile(std::vector<trEntry> & trace, std::string & path);

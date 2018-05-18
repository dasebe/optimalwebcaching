#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <limits>

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

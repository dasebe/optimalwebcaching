#pragma once
#include <vector>
#include <string>
#include <iostream>

// trace entry
struct trEntry {
    uint64_t size;
    double utility;
    uint64_t reqCount;

    trEntry(uint64_t nsize, double util, uint64_t nreqCount)
        : size(nsize),
          utility(util),
          reqCount(nreqCount)
    {
    };

    bool operator <(const trEntry &b)
    {
        return utility < b.utility;
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

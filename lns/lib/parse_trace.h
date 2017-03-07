#include <vector>
#include <string>
#include <lemon/smart_graph.h>

using namespace lemon;

// trace entry
struct trEntry {
    const uint64_t id;
    const uint64_t size;
    const uint64_t origTime;
    double decvar;
    size_t nextSeen;
    bool hasNext;
    int arcId;
    long double utility;

    trEntry(uint64_t nid, uint64_t nsize, uint64_t ntime)
        : id(nid),
          size(nsize),
          origTime(ntime)
    {
        decvar = 0;
        nextSeen = 0;
        hasNext = false;
        arcId = -1;
        utility = 0;
    };
};

// from boost hash combine
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


uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path);

//void createMCF(SmartDigraph & g, std::vector<trEntry> & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int64_t> & supplies);

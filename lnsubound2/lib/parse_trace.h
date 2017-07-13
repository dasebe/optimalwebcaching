#include <vector>
#include <string>
#include <lemon/smart_graph.h>

// uncomment to enable graph debugging:
//#define GDEBUG 1
// uncomment to enable lastSeen debugging:
//#define LDEBUG 1
// uncomment to enable other debugging:
#define ODEBUG 1
//uncomment to enable pi debugging:
//#define PIDEBUG 1
// util for debug
inline void logMessage(std::string m, double x, double y, double z) {
    std::cerr << m << "\t" << x << "\t" << y  << "\t" << z << "\n";
}
#ifdef GDEBUG
#define GLOG(m,x,y,z) logMessage(m,x,y,z)
#else
#define GLOG(m,x,y,z)
#endif
#ifdef LDEBUG
#define LLOG(m,x,y,z) logMessage(m,x,y,z)
#else
#define LLOG(m,x,y,z)
#endif
#ifdef ODEBUG
#define OLOG(m,x,y,z) logMessage(m,x,y,z)
#else
#define OLOG(m,x,y,z)
#endif
#ifdef PIDEBUG
#define PILOG(m,x,y,z) logMessage(m,x,y,z)
#else
#define PILOG(m,x,y,z)
#endif


using namespace lemon;

// trace entry
struct trEntry {
    const uint64_t id;
    const uint64_t size;
    const uint64_t origTime;
    double lcost;
    bool hasNext;
    int64_t outerArcId;
    int64_t innerArcId;
    trEntry(uint64_t nid, uint64_t nsize, uint64_t ntime)
        : id(nid),
          size(nsize),
          origTime(ntime),
          lcost(0),
          hasNext(false),
          outerArcId(-1),
          innerArcId(-1)
    {
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

SmartDigraph::Node createMCF(SmartDigraph & g, std::vector<trEntry> & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int64_t> & supplies);

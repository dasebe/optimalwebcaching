#include <vector>
#include <string>
#include <lemon/smart_graph.h>

using namespace lemon;

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
    double dvar; //for the interval that starts here
    double hit; //for the interval that ends here
    double utility;
    size_t nextSeen;
    int arcId;
    bool hasNext : 1;
    bool active : 1;

    trEntry(uint64_t nid, uint64_t nsize)
        : id(nid),
          size(nsize),
          dvar(0),
          hit(0),
          utility(0),
          nextSeen(0),
          arcId(-1),          
          hasNext(false),
          active(false)
    {
    };
}__attribute__((packed));

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


uint64_t parseTraceFile(std::vector<trEntry> & trace, std::string & path);

inline bool isInEjectSet(const double minUtil, const double maxUtil, const trEntry & curEntry) {
    return( curEntry.utility>=minUtil && curEntry.utility<maxUtil );
}

uint64_t createMCF(SmartDigraph & g, std::vector<trEntry > & trace, uint64_t cacheSize, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<double> & cost, SmartDigraph::NodeMap<int64_t> & supplies, const double minUtil, const double maxUtil);



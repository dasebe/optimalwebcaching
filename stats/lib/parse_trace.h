#include <unordered_set>
#include <string>
#include <unordered_map>
#include <lemon/smart_graph.h>

// trace entry
struct trEntry {
  const uint64_t size;
  size_t nextSeen;
  bool hasParent;
  bool hasChild;
  //  std::unordered_set<size_t> parents;

  trEntry(uint64_t nsize)
    : size(nsize),
    nextSeen(0),
    hasParent(false),
    hasChild(false)
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




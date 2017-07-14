#include <lemon/smart_graph.h>

using namespace lemon;

long double solveMCF(SmartDigraph & g, SmartDigraph::ArcMap<int64_t> & cap, SmartDigraph::ArcMap<long double> & cost, SmartDigraph::NodeMap<int64_t> & supplies, SmartDigraph::ArcMap<int64_t> & flow, int solverPar, SmartDigraph::NodeMap<long double> & lnsPi);

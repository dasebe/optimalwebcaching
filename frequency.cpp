#include <glob.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <map>
#include <math.h>
#include "freq_base.cpp"
using namespace std;

typedef tuple<uint64_t, uint64_t, uint64_t> tuple_t; // objectid, lastr?, size
typedef vector<tuple_t>::iterator vector_iterator_t;


// main
int main (int argc, char* argv[])
{
    // parameters
    if(argc != 5) {
        return 1;
    }
    const char* path = argv[1];
    uint64_t max_size(atoll(argv[2]));
    uint64_t lineCount(atoll(argv[3]));
    uint64_t intervalLen(atoll(argv[4]));

    uint64_t currentLine = 0;

    ifstream infile;
    infile.open(path);

    uint64_t t, id, size;

    StaticCache sc(max_size);
    sc.startStatistics();

    cerr << "running...\n";

    deque<tuple_t> pastDeq;

    while (infile >> t >> id >> size)
    {
        sc.request(id, size);
        pastDeq.push_back(make_tuple(t, id, size));
        if(pastDeq.size() > intervalLen) {
            //            cerr << "antiReq " << currentLine << "\n";
            const tuple_t frontElem = pastDeq.front();
            sc.antiRequest(get<1>(frontElem),get<2>(frontElem));
            pastDeq.pop_front();
        }

        if(++currentLine % intervalLen == 0) {
            cerr << "acHist " << currentLine << "\n";
            sc.activateHistory();
        }
        if(currentLine > lineCount)
            break;
    }
    infile.close();

    cout << "Freq " << max_size << " hitc " << sc.getHits() << " reqc " << currentLine << " OHR " << double(sc.getHits())/currentLine << "\n";
  
    return 0;
}

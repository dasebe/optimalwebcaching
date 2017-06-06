#include <glob.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <math.h>
#include "freq_base.cpp"
using namespace std;

typedef tuple<long, long, long> tuple_t; // objectid, lastr?, size
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

  long currentLine = 0;

  ifstream infile;
  infile.open(path);

  long t, id, size;

  StaticCache sc(max_size);
  sc.startStatistics();

  cerr << "running...\n";

  while (infile >> t >> id >> size)
    {
        sc.request(id, size);
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


#ifndef DATA_H
#define DATA_H

#include "time.h"

#include <vector>
#include <string>

using namespace std;

// structure for holding the data for one timestamp, i.e. one row from the data
// file
struct DataRow {
  int year, month, day, hour, minute, second;
  double B, Bx, By, Bz, Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta;
};

// class for storing, reading and manipulating the data
// generally speaking the amounts of data are huge, it should be taken into
// account here
class Data {
    vector<DataRow> data; // vector of data rows, one row is one timestamp
  public:
    Data(); // constructor
    Data& readFile(string); // read data from file
    Data& filter(int, int); // filter data vector by time (from, until) using
                            // unixtime as limits
    Data& filter(Time, Time); // filter data vector by time (from, until) using
                              // Time objects
    Data copy(); // copy the data object, to be used before filtering
};

#endif


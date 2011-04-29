
#ifndef DATA_H
#define DATA_H

#include "time.h"

#include <eigen3/Eigen/Dense>

#include <vector>
#include <string>

// structure for holding the data for one timestamp, i.e. one row from the
// data file
struct DataRow {
  int year, month, day, hour, minute, second;
  double B, Bx, By, Bz, Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta;
};

// structure for holding the data in the form of Eigen3 vectors
struct DataVectors {
  Eigen::VectorXi year, month, day, hour, minute, second;
  Eigen::VectorXd B, Bx, By, Bz, Vp, Vx, Vy, Vz, Pth, Np, Tp, Vth, beta;
};

// class for storing, reading and manipulating the data
// generally speaking the amounts of data are huge, it should be taken into
// account here
class Data {
    std::vector<DataRow> data; // vector of data rows, one row is one timestamp
    DataVectors vectors; // Eigen3 vectors of data, yes, a copy, but worth it
  public:
    Data& readFile(std::string); // read data from file
    Data& readFile(std::string, Time, Time); // read data file and filter it
    Data& filter(Time, Time); // filter data vector by time (from, until)
                              // using Time objects
    // getters for data, both in rows and columns
    const DataRow& row(int i) const {return data[i];}
    const std::vector<DataRow>& rows() const {return data;}
    const DataVectors& cols() const {return vectors;}
  private:
     // readFile function with optional time limit patameters (NULL pointers)
    void readFile(std::string, Time*, Time*);
    void initVectors(); // initialize Eigen3 vectors of data
};

#endif


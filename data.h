
#ifndef DATA_H
#define DATA_H

// project headers
#include "my_time.h"
#include "axes.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_interp.h>
// standard headers
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
    Data() {}; // default constructor
    Data(DataVectors); // constructor from DataVectors object
    Data& readFile(std::string); // read data from file
    // read data file and filter it
    Data& readFile(std::string, My::Time, My::Time);
    // filter data vector by time (from, until) using Time objects
    Data& filter(My::Time, My::Time);
    // getters for data, both in rows and columns
    const DataRow& row(int i) const {return data[i];}
    const std::vector<DataRow>& rows() const {return data;}
    const DataVectors& cols() const {return vectors;}
    // project data to new coordinates
    Data& project(Axes);
    // resample data
    Data& resample(const int,
                   const gsl_interp_type* interpType = gsl_interp_linear);
  private:
     // readFile function with optional time limit patameters (NULL pointers)
    void readFile(std::string, My::Time*, My::Time*);
    void initVectors(); // initialize Eigen3 vectors of data
};

#endif


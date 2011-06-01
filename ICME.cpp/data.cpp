
#include "data.h"
#include "time.h"
#include "event.h"

#include <eigen3/Eigen/Dense>

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

using namespace std;
using namespace Eigen;

// constructors

// to read data from file and store it in Data object
Data& Data::readFile(string dataFilePath) {
  // read file using NULL pointers for time limits
  readFile(dataFilePath, NULL, NULL);

  return *this; // chained method
}

Data& Data::readFile(string dataFilePath, Time beginTime, Time endTime) {
  // read file using pointers to time limits
  readFile(dataFilePath, &beginTime, &endTime);

  return *this; // chained method
}

// other methods

// slice data using Time objects limits
Data& Data::filter(Time beginTime, Time endTime) {
  int beginIndex, endIndex; // indices for minimum and maximum limits of the
                            // data interval
  Time currentTime; // Time object for storing of the current time

  // iterate through data array
  for (int i=0; i < data.size(); i++) {
    // initialize current Time object
    currentTime = Time(data[i].year, data[i].month, data[i].day,
      data[i].hour, data[i].minute, data[i].second);
    // searching for lower index
    if (currentTime < beginTime) beginIndex = i;
    // searching for upper index
    if (currentTime > endTime) {
      endIndex = i;
      break; // we are already above the upper limit
    }
  } // end of iteration through data

  // erase data outside desired time interval
  data.erase(data.begin()+endIndex, data.end()); // upper part
  data.erase(data.begin(), data.begin()+beginIndex+1); // lower part

  initVectors(); // initialize Eigen3 vectors of data

  return *this; // chained method
}

void Data::readFile(string dataFilePath, Time* beginTime, Time* endTime) {
  ifstream dataFileStream; // stream from the file with data
  string dataFileLine; // a single line from a file as a string
  istringstream dataFileLineStream; // a stream from a line from a file
  DataRow dataRow; // this structure is used to store data for each line from
                   // a data file, i.e. data for a single timestamp
  Time currentTime; // Time object for storing current time for comparison

  // clean data vector, to be on the safe side
  data.clear();
  // open data file as a stream
  dataFileStream.open(dataFilePath.c_str());
  // check if the file was opened correctly
  if (dataFileStream.is_open()) {
    // start iterating through data file line, one timestamp at a time
    while (dataFileStream.good()) {
      // get line from the file as a string
      getline(dataFileStream, dataFileLine);
      // check if the line contains actually data
      if (!dataFileLine.empty() && // check that it's not empty
          dataFileLine[0] != '#') // and not commented out
      {
        // clear string stream from the previus line of data
        dataFileLineStream.clear();
        // use the next line of data as a source for the string stream
        dataFileLineStream.str(dataFileLine);
        // parse the data from the stream of line of the data file
        dataFileLineStream >> dataRow.year >> dataRow.month >> dataRow.day >>
          dataRow.hour >> dataRow.minute >> dataRow.second >> dataRow.B >>
          dataRow.Bx >> dataRow.By >> dataRow.Bz >> dataRow.Vp >>
          dataRow.Vx >> dataRow.Vy >> dataRow.Vz >> dataRow.Pth >>
          dataRow.Np >> dataRow.Tp >> dataRow.Vth >> dataRow.beta;
        // do filtering only if pointers to time limits are not NULL
        if (beginTime != NULL && endTime != NULL) {
          // initialize current Time object with time data
          currentTime = Time(dataRow.year, dataRow.month, dataRow.day,
            dataRow.hour, dataRow.minute, dataRow.second);
          if (currentTime < *beginTime) { // before the minimum time limit
            continue; // miss it
          } else {
            // converting data to SI units
            dataRow.B   = dataRow.B*1e-9;  // nT    -> T
            dataRow.Bx  = dataRow.Bx*1e-9; // nT    -> T
            dataRow.By  = dataRow.By*1e-9; // nT    -> T
            dataRow.Bz  = dataRow.Bz*1e-9; // nT    -> T
            dataRow.Vp  = dataRow.Vp*1e3;  // km/s  -> m/s
            dataRow.Vx  = dataRow.Vx*1e3;  // km/s  -> m/s
            dataRow.Vy  = dataRow.Vy*1e3;  // km/s  -> m/s
            dataRow.Vz  = dataRow.Vz*1e3;  // km/s  -> m/s
            dataRow.Pth = dataRow.Pth*1e-9; // nPa   -> Pa
            dataRow.Np  = dataRow.Np*1e6;  // 1/cm3 -> 1/m3
            dataRow.Vth = dataRow.Vth*1e3; // km/s  -> m/s
            data.push_back(dataRow); // push the row of data to the data vector
          }
          if (currentTime > *endTime) break; // after the maximum time limit
        }
      }
    } // end of iteration through the lines of the data file
    initVectors(); // initialize Eigen3 vectors of data
  }
}

// project data to new coordinates
Data& Data::project(Axes axes) {

  Vector3d V, B; // initialize temporary vectors for speed and magnetic field

  // iterate through data and project it to new axes
  for (int i=0; i < data.size(); i++) {
    // fill temporary vectors
    V(0) = data[i].Vx;
    V(1) = data[i].Vy;
    V(2) = data[i].Vz;
    B(0) = data[i].Bx;
    B(1) = data[i].By;
    B(2) = data[i].Bz;
    // project data to new axes
    data[i].Vx = V.dot(axes.x);
    data[i].Vy = V.dot(axes.y);
    data[i].Vz = V.dot(axes.z);
    data[i].Bx = B.dot(axes.x);
    data[i].By = B.dot(axes.y);
    data[i].Bz = B.dot(axes.z);
  } // end of iteration through the data

  initVectors(); // recalculate Eigen3 vectors of data

  return *this; // chained method
}

// initialize Eigen3 vectors of data. very stupid function... need to find
// a way to get rid of this variable lists
void Data::initVectors() {
  int n = data.size(); // length of data vectors
  // resize dynamic vectors for data
  vectors.year.resize(n);
  vectors.month.resize(n);
  vectors.day.resize(n);
  vectors.hour.resize(n);
  vectors.minute.resize(n);
  vectors.second.resize(n);
  vectors.B.resize(n);
  vectors.Bx.resize(n);
  vectors.By.resize(n);
  vectors.Bz.resize(n);
  vectors.Vp.resize(n);
  vectors.Vx.resize(n);
  vectors.Vy.resize(n);
  vectors.Vz.resize(n);
  vectors.Pth.resize(n);
  vectors.Np.resize(n);
  vectors.Tp.resize(n);
  vectors.Vth.resize(n);
  vectors.beta.resize(n);
  // fill the vectors with data
  for (int i=0; i < n; i++) {
    vectors.year(i) = data[i].year;
    vectors.month(i) = data[i].month;
    vectors.day(i) = data[i].day;
    vectors.hour(i) = data[i].hour;
    vectors.minute(i) = data[i].minute;
    vectors.second(i) = data[i].second;
    vectors.B(i) = data[i].B;
    vectors.Bx(i) = data[i].Bx;
    vectors.By(i) = data[i].By;
    vectors.Bz(i) = data[i].Bz;
    vectors.Vp(i) = data[i].Vp;
    vectors.Vx(i) = data[i].Vx;
    vectors.Vy(i) = data[i].Vy;
    vectors.Vz(i) = data[i].Vz;
    vectors.Pth(i) = data[i].Pth;
    vectors.Np(i) = data[i].Np;
    vectors.Tp(i) = data[i].Tp;
    vectors.Vth(i) = data[i].Vth;
    vectors.beta(i) = data[i].beta;
  }
}


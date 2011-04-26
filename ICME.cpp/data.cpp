
#include "data.h"
#include "time.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

using namespace std;

// constructors

// TODO: resolve duplication of code in constructors!

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
  data.erase(data.begin(), data.begin()+beginIndex); // lower part
  data.erase(data.begin()+endIndex, data.end()); // upper part

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
            data.push_back(dataRow); // push the row of data to the data vector
          }
          if (currentTime > *endTime) break; // after the maximum time limit
        }
      }
    } // end of iteration through the lines of the data file
  }
}


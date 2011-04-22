
#include "data.h"
#include "time.h"

#include <fstream>
#include <sstream>
#include <string>

using namespace std;

// constructor
Data::Data() {}

// to read data from file and store it in Data object
Data& Data::readFile(string dataFilePath) {
  ifstream dataFileStream; // stream from the file with data
  string dataFileLine; // a single line from a file as a string
  istringstream dataFileLineStream; // a stream from a line from a file
  DataRow dataRow; // this structure is used to store data for each line from a
                   // data file, i.e. data for a single timestamp

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
          dataRow.Bx >> dataRow.By >> dataRow.Bz >> dataRow.Vp >> dataRow.Vx >>
          dataRow.Vy >> dataRow.Vz >> dataRow.Pth >> dataRow.Np >> dataRow.Tp >>
          dataRow.Vth >> dataRow.beta;
        // push the row of data to the data vector
        data.push_back(dataRow);
      }
    } // end of iteration through the lines of the data file
  }

  return *this; // chained method
}

// slice data using unixtime limits
Data& filter(int beginUnixtime, int endUnixtime) {

}

// slice data using Time objects limits
Data& filter(Time beginTime, Time endTime) {

}


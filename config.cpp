
// project headers
#include "config.h"
#include "time.h"
// standard headers
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <iostream>

using namespace std;
using namespace My;

// constructor
Config::Config() {}

// read config data from file
Config& Config::readFile(string configFilePath) {
  ifstream configFileStream; // file stream
  string configFileLine; // file line, a string
  istringstream configFileLineStream; // string stream for a file line
  ConfigRow configRow; // structure for storing the data from a single line
                       // of a file
  string beginDate, endDate, beginTime, endTime; // begin and end time strings

  // clear the data vector, just in case
  data.clear();
  // open file as an input stream
  configFileStream.open(configFilePath.c_str());
  //check if everything is ok
  if (configFileStream.is_open()) {
    // start iterating through file, quit when encounter end of the file
    while (configFileStream.good()) {
      // get a line from the file as a string
      getline(configFileStream, configFileLine);
      // we continue only if the line is not empty and not a comment
      if (!configFileLine.empty() && // check if the line is empty
          configFileLine[0] != '#' && // check if it's commented
          configFileLine[0] != '!' && // check if it's the first line
          configFileLine[0] != ' ')   // check if it's an empty line
      {
        // treating each line as a string stream
        // clean the string stream, required on every iteration
        configFileLineStream.clear();
        // set the string from the file as the source for the stream
        configFileLineStream.str(configFileLine);
        // get the columned data from the string stream and write it to the
        // ConfigRow structure
        configFileLineStream >> configRow.flag >> configRow.quality >>
          configRow.toGsr >> configRow.toMva >> configRow.toCm >>
          configRow.toTm >> configRow.toHm >> configRow.toOm >>
          configRow.toSave >> configRow.spacecraft >> beginDate >> beginTime >>
          endDate >> endTime >> configRow.samplingInterval >>
          configRow.fittingFuntion >> configRow.fittingParameterBdr >>
          configRow.fittingParameterCtr >> configRow.Nx >> configRow.ratio >>
          configRow.minY >> configRow.maxY >> configRow.order;
        // calculate begin and end Time objects instead of strings
        configRow.beginTime = Time(beginDate+' '+beginTime);
        configRow.endTime = Time(endDate+' '+endTime);
        // push the ConfigRow structure to the vector of config data
        data.push_back(configRow);
      }
    } // end of iteration through the config file
  }

  return *this; // chained method
}

// filter config data by quality and by execution flag
// we are filtering in place, i.e. by removing events which doesn't pass
// filtering, not by copying an object
Config& Config::filter(char eventQuality) {
  vector<ConfigRow>::iterator configIter; // iterator for going through the
                                          // vector of config data
  // starting the iteration
  for (configIter=data.begin(); configIter != data.end(); configIter++) {
    // check if the event is not of sufficient quality or is not intended for
    // analysis
    if (configIter->quality != eventQuality || configIter->flag == 0) {
      configIter = data.erase(configIter); // erase the event from the data
      configIter--; // after erasing we get the iterator pointing to the next
                    // element, but we don't need it, so we make one step back
    }
  } // end of the iteration

  return *this; // chained method
}


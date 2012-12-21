
#ifndef CONFIG_H
#define CONFIG_H

// project headers
#include "my_time.h"
// standard headers
#include <string>
#include <vector>

// this structure represents a single row from the config file
struct ConfigRow {
  int flag, Nx, order;
  double samplingInterval;
  char quality;
  bool toGsr, toMva, toCm, toTm, toHm, toOm, toSave;
  std::string spacecraft, fittingFuntion;
  double ratio, minY, maxY, fittingParameterCtr, fittingParameterBdr,
         minTheta, maxTheta, minPhi, maxPhi;
  My::Time beginTime, endTime;
};

// config class, reads the data from file and filters it
class Config {
    // here the array of config file rows is stored
    std::vector<ConfigRow> data;
  public:
    Config(); // constructor
    Config& readFile(std::string); // read config data from file
    Config& filter(char); // filter the config data
    // use this to get one row of the config data, i.e. config for one event
    ConfigRow row(int i) const {return data[i];}
    // use this to get all the rows of the config data, i.e. for all events
    std::vector<ConfigRow> rows() const {return data;}
};

#endif


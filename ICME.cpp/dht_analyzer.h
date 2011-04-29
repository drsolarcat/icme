
#ifndef DHT_ANALYZER_H
#define DHT_ANALYZER_H

#include "event.h"

#include <eigen3/Eigen/Dense>

// this class is used for searching for the deHoffmann-Teller frame
class DhtAnalyzer {
  public:
    void analyze(Event&); // main trigger for dHT analysis
  private:
    // loop through the range of speeds and choose one most close to the
    // deHoffmann-Teller speed
    Eigen::Vector3d loop(Event&, double, double, double,
                                 double, double, double,
                                 double, double, double);
};

#endif


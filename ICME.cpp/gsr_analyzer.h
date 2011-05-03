
#ifndef GSR_ANALYZER_H
#define GSR_ANALYZER_H

#include "event.h"

#include <eigen3/Eigen/Dense>

// this class is used to carry GSR analysis
class GsrAnalyzer {
  public:
    void analyze(Event&); // main trigger to carry the analysis,
                          // the results are stored in the Event object
  private:
    // loop through angles and find the optimal angles for new axes
    GsrRun loopAxes(Event&, double, double, double,
                            double, double, double);

};

#endif


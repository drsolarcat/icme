
#ifndef GSR_ANALYZER_H
#define GSR_ANALYZER_H

// project headers
#include "event.h"
#include "my_time.h"
// library headers
#include <eigen3/Eigen/Dense>
// standard headers
#include <string>

// this class is used to carry GSR analysis
class GsrAnalyzer {

  public:
    void analyze(Event&); // main trigger to carry the analysis,
                          // the results are stored in the Event object
  private:
    // loop through angles and find the optimal angles for new axes
    GsrResults loopAxes(Event&, double, double, double,
                                double, double, double);
    // detect axes if the flux rope
    GsrResults& detectAxes(Event&, GsrResults&);
    // detect flux rope axes in the non-symmetrical case
    GsrResults& detectAxesByResidue(Event&, GsrResults&);
    // detect flux rope axes in case of symmetry
    GsrResults& detectAxesByVariance(Event&, GsrResults&);
    // compute magnetic field map for a given GSR run
    GsrResults& computeMap(Event&, GsrResults&);
    // get transformation matrix
    Eigen::Matrix3d getTransformationMatrix(Event&, std::string, My::Time);
};

#endif


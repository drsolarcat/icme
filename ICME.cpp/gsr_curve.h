
#ifndef GSR_CURVE_H
#define GSR_CURVE_H
// project headers
#include "branched_curve.h"
#include "event.h"
#include "axes.h"
// library headers
#include <eigen3/Eigen/Dense>

// this class of branched curves supports initializing by vector potential and
// transverse pressure, specific to GSR method
class GsrCurve: public BranchedCurve {
  public:
    GsrCurve();
    GsrCurve(Event&, Axes); // physical constructor
};

#endif


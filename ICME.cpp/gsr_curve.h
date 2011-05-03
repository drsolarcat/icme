
#ifndef GSR_CURVE_H
#define GSR_CURVE_H

#include "branched_curve.h"
#include "event.h"
#include "axes.h"

#include <eigen3/Eigen/Dense>

// this class of branched curves supports initializing by vector potential and
// transverse pressure, specific to GSR method
class GsrCurve: public BranchedCurve {
  public:
    GsrCurve(Event&, Axes); // physical constructor
};

#endif



#ifndef BRANCHED_CURVE_H
#define BRANCHED_CURVE_H

#include "curve.h"

#include <vector>

// this class of curves supports branching
class BranchedCurve: public Curve {
  protected:
    bool c_isBranched; // weather it is branched
    std::vector<Curve> c_branches; // holds to branch curve, if any
  public:
    std::vector<Curve>& branches() {return c_branches;} // branches accessor
    bool isBranched() {return c_isBranched;} // true if branched
    std::vector<Curve>& initBranches(); // determine the branches of the curve
};

#endif


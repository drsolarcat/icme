
#ifndef BRANCHED_CURVE_H
#define BRANCHED_CURVE_H

#include "curve.h"

#include <vector>

// this class of curves supports branching
class BranchedCurve: public Curve {
  protected:
    bool _isBranched; // whether it is branched
    int _minLeftIndex, _maxIndex, _minRightIndex;
    double _residue;
    double _combinedResidue;
    int _branchLength;
    std::vector<Curve> _branches; // holds to branch curve, if any
  public:
    // accessors for boundary indices
    const int minLeftIndex() const {return _minLeftIndex;}
    const int maxIndex() const {return _maxIndex;}
    const int minRightIndex() const {return _minRightIndex;}
    std::vector<Curve>& branches() {return _branches;} // branches accessor
    bool isBranched() {return _isBranched;} // true if branched
    double residue() {return _residue;} // residue accessor
    double combinedResidue() {return _combinedResidue;} // residue accessor
    double branchLength() {return _branchLength;} // residue accessor
    BranchedCurve& initBranchesByTracing();
    BranchedCurve& initBranchesByExtremums();
    BranchedCurve& initBranches(); // determine the branches of the curve
    BranchedCurve& initBranches(int, int, int); // initialize curves by
                                                // boundary indices
    BranchedCurve& computeResidue(); // compute the residue between
                                     // the branches of the curve
};

#endif


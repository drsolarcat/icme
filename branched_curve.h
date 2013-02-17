
#ifndef BRANCHED_CURVE_H
#define BRANCHED_CURVE_H

// project headers
#include "curve.h"
// standard headers
#include <vector>
#include <string>

// this class of curves supports branching
class BranchedCurve: public Curve {
  protected:
    bool _isBranched; // whether it is branched
    // boudary indices for the branches
    int _leftIndex, _centerIndex, _rightIndex;
    double _originalResidue; // original residue as in Hu, Sonnerup (1999)
    double _combinedResidue; // combined residue as in Isavnin et al. (2011)
    int _branchLength; // length of a shorter branch
    std::vector<Curve> _branches; // holds to branch curve, if any
  public:
    BranchedCurve(Eigen::VectorXd x, Eigen::VectorXd y) {
      _vectors.x = x;
      _vectors.y = y;
    }
    BranchedCurve() {}
    // accessors for boundary indices
    const int leftIndex() const {return _leftIndex;}
    const int centerIndex() const {return _centerIndex;}
    const int rightIndex() const {return _rightIndex;}
    // branches accessor
    const std::vector<Curve>& branches() const {return _branches;}
    const bool isBranched() const {return _isBranched;} // true if branched
    // residue accessors
    const double originalResidue() const {return _originalResidue;}
    const double combinedResidue() const {return _combinedResidue;}
    // branch length accessor
    const double branchLength() const {return _branchLength;}
    // init branches, general method
    BranchedCurve& initBranches(std::string);
    // init branches by boundary indices
    BranchedCurve& initBranches(int, int, int);
    // compute resudual parameters
    BranchedCurve& computeResidue();
  protected:
    // init branches by tracing
    BranchedCurve& initBranchesByTracing();
    // init branches by extremum analysis
    BranchedCurve& initBranchesByExtremums();
};

#endif


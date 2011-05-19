
#include "branched_curve.h"

#include <eigen3/Eigen/Dense>

#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

// initialize branches of the curve
BranchedCurve& BranchedCurve::initBranches() {

  // copy curve data vectors
  VectorXd x = _vectors.x;
  VectorXd y = _vectors.y;

  // smooth the curve data
  filterRunningAverage(x, 5);
  filterRunningAverage(y, 5);

  // let X be always with positive camber
  if (abs(x.minCoeff()) > abs(x.maxCoeff())) x = -x;

  // point of X maximum
  int maxIndexX;
  x.maxCoeff(&maxIndexX);

  // initial guess for branches limits
  _isBranched = true;
  _minLeftIndex = 0;
  _maxIndex = maxIndexX;
  _minRightIndex = x.size();

  // if X maximum is at the boundaries of X then no branches
  if (maxIndexX == 0 || maxIndexX == x.size()-1) {
    _isBranched = false;
  } else {
    // initialize minimums on both sides of X maximum
    int minLeftIndexX, minRightIndexX;
    double minLeftX = x.head(maxIndexX+1).minCoeff(&minLeftIndexX);
    double minRightX = x.tail(x.size()-maxIndexX).minCoeff(&minRightIndexX);
    minRightIndexX += maxIndexX;

    // align the minimums, so that they are equal Y-valued
    // if the left minimum is higher than the right one
    if (minLeftX >= minRightX)) {
      // go from the maximum to the right minimum
      for (int i = maxIndexX; i < x.size(); i++) {
        if (x(i) < minLeftX) {
          minRightIndexX = i-1;
          break;
        }
      }
    } else { // if the right minimum is higher than the left one
      // go from the maximum to the left minimum
      for (int i = maxIndexX; i >= 0; i--) {
        if (x(i) < minRightX) {
          minLeftIndexX = i+1;
          break;
        }
      }
    }

    // save branches boundaries
    _minLeftIndex = minLeftIndexX;
    _minRightIndex = minRightIndexX;

//    cout << _minLeftIndex << ':' << _maxIndex << ':' << _minRightIndex << endl;
  }

  // initialize branch curves
  if (_isBranched) {
    initBranches(_minLeftIndex, _maxIndex, _minRightIndex);
  }

  return *this; // chained method
}

// initialize branches by the boudaries
BranchedCurve& BranchedCurve::initBranches(int minLeftIndex,
                                           int maxIndex,
                                           int minRightIndex)
{
  _branches.clear(); // clear branches
  // initialize left minimum to maximum branch
  _branches.push_back(
    Curve(_vectors.x.segment(minLeftIndex, maxIndex-minLeftIndex),
          _vectors.y.segment(minLeftIndex, maxIndex-minLeftIndex)));
  // initialize maximum to right minimum branch
  _branches.push_back(
    Curve(_vectors.x.segment(maxIndex, minRightIndex-maxIndex),
          _vectors.y.segment(maxIndex, minRightIndex-maxIndex)));

  return *this; // chained method
}

// compute residues of the curve branches
BranchedCurve& BranchedCurve::computeResidue() {

}


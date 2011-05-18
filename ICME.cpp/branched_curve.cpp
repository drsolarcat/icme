
#include "branched_curve.h"

#include <eigen3/Eigen/Dense>

#include <iostream>
#include <algorithm>

using namespace std;
using namespace Eigen;

// initialize branches of the curve
BranchedCurve& BranchedCurve::initBranches() {

  VectorXd x = _vectors.x;
  VectorXd y = _vectors.y;

  filterRunningAverage(x, 5);
  filterRunningAverage(y, 5);

  if (abs(x.minCoeff()) > abs(x.maxCoeff())) x = -x;

  int maxIndexX;
  x.maxCoeff(&maxIndexX);

  _isBranched = true;
  _minLeftIndex = 0;
  _maxIndex = maxIndexX;
  _minRightIndex = x.size();

  if (maxIndexX == 0 || maxIndexX == x.size()-1) {
    _isBranched = false;
  } else {
    int minLeftIndexX, minRightIndexX;
    double minLeftX = x.head(maxIndexX+1).minCoeff(&minLeftIndexX);
    double minRightX = x.tail(x.size()-maxIndexX).minCoeff(&minRightIndexX);
    minRightIndexX += maxIndexX;

    if (minLeftX == max(minLeftX, minRightX)) {
      for (int i = maxIndexX; i < x.size(); i++) {
        if (x(i) < minLeftX) {
          minRightIndexX = i-1;
          break;
        }
      }
    } else {
      for (int i = maxIndexX; i >= 0; i--) {
        if (x(i) < minRightX) {
          minLeftIndexX = i+1;
          break;
        }
      }
    }

    _minLeftIndex = minLeftIndexX;
    _minRightIndex = minRightIndexX;

    cout << _minLeftIndex << ':' << _maxIndex << ':' << _minRightIndex << endl;
  }

  if (_isBranched) {
    initBranches(_minLeftIndex, _maxIndex, _minRightIndex);
  }

  return *this; // chained method
}

BranchedCurve& BranchedCurve::initBranches(int minLeftIndex,
                                           int maxIndex,
                                           int minRightIndex)
{
  _branches.clear();
  _branches.push_back(
    Curve(_vectors.x.segment(minLeftIndex, maxIndex-minLeftIndex),
          _vectors.y.segment(minLeftIndex, maxIndex-minLeftIndex)));
  _branches.push_back(
    Curve(_vectors.x.segment(maxIndex, minRightIndex-maxIndex),
          _vectors.y.segment(maxIndex, minRightIndex-maxIndex)));
}

BranchedCurve& BranchedCurve::computeResidue() {

}


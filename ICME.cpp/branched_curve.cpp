
// project headers
#include "branched_curve.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_version.h>
// standard headers
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;
using namespace Eigen;

// init branches, general method
BranchedCurve& BranchedCurve::initBranches(string method) {
  if (method == "tracing") {
    initBranchesByTracing(); // init branches by tracing
  } else if (method == "extremums") {
    initBranchesByExtremums(); // init branches by extremums
  }
}

// init branches by tracing
BranchedCurve& BranchedCurve::initBranchesByTracing() {

  // copy the curve data vectors
  VectorXd x = _vectors.x;
  VectorXd y = _vectors.y;

  // smooth the curve, we are smoothing only the Y-part (Pt in GSR)
  filterRunningAverage(y, 0.07*size());
  filterRunningAverage(y, 0.1*size());
//  filterRunningAverage(y, 5);

  // initial values for boundary indices, stored locally
  int leftIndex = 0,
      centerIndex,
      rightIndex = size()-1;
  x.maxCoeff(&centerIndex);

  // initial values for boundary indices, stored in a curve object
  _centerIndex = centerIndex;
  _leftIndex = leftIndex;
  _rightIndex = rightIndex;
  _isBranched = false;

  // mirror the X if the shape is not concave
  if (x(centerIndex) < x(leftIndex) &&
      x(centerIndex) < x(rightIndex))
  {
    x = -x; // mirror X
  }

  // criteron for coincident point of the branches
  const double eps2 = pow(0.0001*(x.maxCoeff()-x.minCoeff()), 2)+
                      pow(0.01*(y.maxCoeff()-y.minCoeff()), 2); // TODO: too big

  // if the central index is not on the edge
  if (x(centerIndex) != x(rightIndex) &&
      x(centerIndex) != x(leftIndex))
  {
//    while (rightIndex-leftIndex > 0.4*size()) {
    // search until left and right indices coincide
    while (rightIndex > leftIndex) {
      // move to the central indice step by step
      if (x(leftIndex) < x(rightIndex)) {
        leftIndex++; // step from left to right
      } else if (x(leftIndex) > x(rightIndex)) {
        rightIndex--; // step from right to left
      }
      // check the criterion: how close are the points corresponding
      // to the boundaries are?
      if (pow(x(leftIndex)-x(rightIndex), 2)+
          pow(y(leftIndex)-y(rightIndex), 2) < eps2)
      {
        // close enough, save indices and break
        _rightIndex = rightIndex;
        _leftIndex = leftIndex;
        _centerIndex = centerIndex;
        _isBranched = true;
        break;
      }
    }
  }

  // initialize branch curves
  if (_isBranched) {
    initBranches(_leftIndex, _centerIndex, _rightIndex);
  }

//  cout << _minLeftIndex << ':' << _maxIndex << ':' << _minRightIndex << endl;

  return *this; // chained method
}

// initialize branches of the curve
BranchedCurve& BranchedCurve::initBranchesByExtremums() {

  // copy curve data vectors
  VectorXd x = _vectors.x;
  VectorXd y = _vectors.y;

  // smooth the curve data
  filterRunningAverage(y, 0.07*size());
  filterRunningAverage(y, 0.1*size());
//  filterRunningAverage(y, 5);

  // let X be always with positive camber
  if (abs(x.minCoeff()) > abs(x.maxCoeff())) x = -x;

  // point of X maximum
  int maxIndexX;
  x.maxCoeff(&maxIndexX);

  // initial guess for branches limits
  _isBranched = true;
  _leftIndex = 0;
  _centerIndex = maxIndexX;
  _rightIndex = x.size();

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
    if (minLeftX >= minRightX) {
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

    const int delta = int(0.05*size());

    // save branches boundaries
    _leftIndex = minLeftIndexX;
    _rightIndex = minRightIndexX;
  }

  // initialize branch curves
  if (_isBranched) {
    initBranches(_leftIndex, _centerIndex, _rightIndex);
  }

//  cout << _minLeftIndex << ':' << _maxIndex << ':' << _minRightIndex << endl;

  return *this; // chained method
}

// initialize branches by the boudaries
BranchedCurve& BranchedCurve::initBranches(int leftIndex,
                                           int centerIndex,
                                           int rightIndex)
{
  _branches.clear(); // clear branches
  // initialize left minimum to maximum branch
  _branches.push_back(
    Curve(_vectors.x.segment(leftIndex, centerIndex-leftIndex+1),
          _vectors.y.segment(leftIndex, centerIndex-leftIndex+1)));

  // initialize maximum to right minimum branch
  _branches.push_back(
    Curve(_vectors.x.segment(centerIndex, rightIndex-centerIndex+1),
          _vectors.y.segment(centerIndex, rightIndex-centerIndex+1)));

  return *this; // chained method
}

// compute residues of the curve branches
BranchedCurve& BranchedCurve::computeResidue() {

  // minimal length of the branch
  _branchLength = min(_centerIndex - _leftIndex + 1,
                      _rightIndex - _centerIndex + 1);

  // interpolation type used for resampling
  const gsl_interp_type* interp_type = gsl_interp_linear;

  /* if GSL_VERSION >= 1.15 */
  /* if (_isBranched && m >= gsl_interp_type_min_size(interp_type)) { */
  /* else */
  const int interp_min_size = 2;
  if (_isBranched && _branchLength >= interp_min_size) {
  /* endif */
    // copy the curve branches
    Curve curveIn(_branches[0]);
    Curve curveOut(_branches[1]);

    // init the limits for the resampled branches
    double minX = max(curveIn.cols().x.minCoeff(),
                      curveOut.cols().x.minCoeff());
    double maxX = min(curveIn.cols().x.maxCoeff(),
                      curveOut.cols().x.maxCoeff());

    // resample the branches
    curveIn.resample(minX, maxX, 2*size(), interp_type);
    curveOut.resample(minX, maxX, 2*size(), interp_type);

    // init the limits of the Y values of the curves
    double minY = min(curveIn.cols().y.minCoeff(),
                      curveOut.cols().y.minCoeff());
    double maxY = max(curveIn.cols().y.maxCoeff(),
                      curveOut.cols().y.maxCoeff());

    // compute residue
    // original residue as in Hu, Sonnerup (1999)
    _originalResidue = (curveIn.cols().y-curveOut.cols().y).norm()/abs(maxY-minY);
    // combined residue as in Isavnin et al. (2011)
    _combinedResidue = _originalResidue*size()/2/_branchLength;
  } else { // if the curve is not branched
    _originalResidue = numeric_limits<double>::infinity();
    _combinedResidue = numeric_limits<double>::infinity();
  }
}



#include "curve.h"

#include <eigen3/Eigen/Dense>

using namespace Eigen;

// default constructor for the Curve objects
Curve::Curve() {}

// smooth curve component by running average filter
Curve& Curve::filterRunningAverage(int span, char axis) {
  if (axis == 'x') { // filter x components
    filterRunningAverage(c_vectors.x, span);
  } else if (axis == 'y') { // filter y component
    filterRunningAverage(c_vectors.y, span);
  }

  return *this; // chained method
}

// return smoothed by running average version of the data vector
// static member function
VectorXd Curve::filteredRunningAverage(VectorXd y, int span) {
  const int delta = (span-1)/2; // number of points to the left and right of
                               // the current point of the data being filtered
  int d = 0; // the same, will be recalculated on every step
  const int n = y.size(); // size of the data vector to be filtered
  VectorXd yy = VectorXd::Zero(n); // new vector, will hold the fitered data

  // begin iteration through the data
  for (int i = 0; i < n; i++) {
    if (i < delta) { // cose to the beginning of the data vector
      d = i;
    } else if (n-i < delta) { // close to the end of the data vector
      d = n - i;
    } else { // somewhere in the middle of the data vector
      d = delta;
    }
    yy(i) = y.segment(i-d,2*d).sum()/(2*d+1); // perform filtering
  } // end iteration through the data

  return yy; // return new filtered data vector
}

// perform filtering by running average inplace
// static member function
void Curve::filterRunningAverage(VectorXd& y, int span) {
  y = Curve::filteredRunningAverage(y, span); // get filtered data vector and
                                              // assign it to the old variable
}


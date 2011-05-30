
#include "curve.h"

#include <eigen3/Eigen/Dense>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>

using namespace Eigen;

// default constructor for the Curve objects
Curve::Curve() {}

// construct the curve out of data vectors
Curve::Curve(VectorXd x, VectorXd y) {
  _vectors.x = x;
  _vectors.y = y;
}

// smooth curve component by running average filter
Curve& Curve::filterRunningAverage(int span, char axis) {
  if (axis == 'x') { // filter x components
    filterRunningAverage(_vectors.x, span);
  } else if (axis == 'y') { // filter y component
    filterRunningAverage(_vectors.y, span);
  }

  return *this; // chained method
}

// resample curve using min and max X limits and number of points
Curve& Curve::resample(double minX, double maxX, const int m) {
  const int n = _vectors.x.size();

//  if (n >= gsl_interp_type_min_size(gsl_interp_linear)) {
  if (n >= 2) {
    // initialize the accelerator object
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    // initialize the interpolating spline
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_linear, n);

    double x[n], y[n]; // arrays for storing vector data
    size_t p[n]; // permuations arrays

    // write curve data to arrays
    for (int i = 0; i < size(); i++) {
      x[i] = _vectors.x(i);
      y[i] = _vectors.y(i);
    }
    gsl_sort_index(p, x, 1, n); // sort curve data by ASC X, get permutations
    gsl_permute(p, x, 1, n); // apply permutations
    gsl_permute(p, y, 1, n); // apply permutations
    gsl_spline_init(spline, x, y, n); // initialize the interpolating spline
    double dx = (maxX-minX)/(m-1); // compute the X step
    _vectors.x.resize(m); // resize curve vector for the resampled curve
    _vectors.y.resize(m); // resize curve vector for the resampled curve
    // save new resampled curve data
    for (int i = 0; i < m; i++) {
       // resampled X
      _vectors.x(i) = minX+i*dx;
      // resampled Y
      _vectors.y(i) = gsl_spline_eval(spline, _vectors.x(i), acc);
    }
    gsl_spline_free (spline); // free spline memory
    gsl_interp_accel_free (acc); // free accelerator memory
  }
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


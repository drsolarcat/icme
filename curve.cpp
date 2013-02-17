
// project header
#include "curve.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
// standard headers
#include <iostream>

using namespace Eigen;
using namespace std;

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
Curve& Curve::resample(const double minX, const double maxX, const int m,
                       const gsl_interp_type* interpType)
{
  sort(); // sort the curve
  unique(); // and remove duplicates

  const int n = _vectors.x.size(); // curve length

  // check if there are enough points for interpolation
  if (n >= gsl_interp_type_min_size(interpType)) {
    // initialize the accelerator object
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    // initialize the interpolating spline
    gsl_spline *spline = gsl_spline_alloc(interpType, n);

    double *x, *y; // arrays for storing vector data

    // write curve data to arrays
    x = _vectors.x.data();
    y = _vectors.y.data();

    gsl_spline_init(spline, x, y, n); // initialize the interpolating spline
    double dx = (maxX-minX)/(m-1); // compute the X step
    _vectors.x.resize(m); // resize curve vector for the resampled curve
    _vectors.y.resize(m); // resize curve vector for the resampled curve
    // save new resampled curve data
    for (int i = 0; i < m; i++) {
      // resampled X
      _vectors.x(i) = (i == m-1 ? maxX : minX+i*dx); // sometimes this tends to
                                                     // go > maxX because of
                                                     // dx precision, it's a fix
      // resampled Y
      _vectors.y(i) = gsl_spline_eval(spline, _vectors.x(i), acc);
    }
    gsl_spline_free(spline); // free spline memory
    gsl_interp_accel_free(acc); // free accelerator memory
  }

  return *this; // chained method
}

// resample the curve using the number of points only
Curve& Curve::resample(const int m,
                       const gsl_interp_type* interpType)
{
  // resampling using X boundaries
  resample(_vectors.x.minCoeff(), _vectors.x.maxCoeff(), m, interpType);

  return *this; // chained method
}

// sort the curve by X data
Curve& Curve::sort() {
  // length of the curve
  const int n = _vectors.x.size();

  double *xArr, *yArr; // arrays for curve data
  size_t p[n]; // permutations array

  // initialize arrays with vectors data
  // WARNING: operating with pointers
  xArr = _vectors.x.data();
  yArr = _vectors.y.data();

  // get permutations array for sorting the X data
  gsl_sort_index(p, xArr, 1, n);
  // sort the data
  gsl_permute(p, xArr, 1, n);
  gsl_permute(p, yArr, 1, n);

  // initialize new sorted X and Y vectors
  Map<VectorXd> xSorted(xArr, n);
  Map<VectorXd> ySorted(yArr, n);

  // substitute old data vectors with new ones
  _vectors.x = xSorted;
  _vectors.y = ySorted;

  return *this; // chained method
}

// remove X duplicates to the curve into a true function
Curve& Curve::unique() {
  // length of the curve
  const int n = _vectors.x.size();

  // initialize temporary vector objects to store unique data
  VectorXd xUnique = VectorXd::Zero(n),
           yUnique = VectorXd::Zero(n);

  int k = 0; // new index initialization
  for (int i = 0; i < n; i++) { // begin iteration through data
    if (i >= 1 && _vectors.x(i) == _vectors.x(i-1)) { // if not unique
      k--; // step back
      xUnique(k) = _vectors.x(i); // save X
      yUnique(k) = (_vectors.y(i)+_vectors.y(i-1))/2; // average Y
    } else { // if unique
      xUnique(k) = _vectors.x(i);
      yUnique(k) = _vectors.y(i);
    }
    k++; // step forward
  }

  // resize vectors
  xUnique.conservativeResize(k);
  yUnique.conservativeResize(k);

  // substitute old data vectors with unique ones
  _vectors.x = xUnique;
  _vectors.y = yUnique;

  return *this; // chained method
}

// return smoothed by running average version of the data vector
// static member function
VectorXd Curve::filteredRunningAverage(const VectorXd& y, const int span) {
  const int delta = (span-1)/2; // number of points to the left and right of
                               // the current point of the data being filtered
  int d = 0; // the same, will be recalculated on every step
  const int n = y.size(); // size of the data vector to be filtered
  VectorXd yy = VectorXd::Zero(n); // new vector, will hold the fitered data

  // begin iteration through the data
  for (int i = 0; i < n; i++) {
    if (i < delta) { // close to the beginning of the data vector
      d = i;
    } else if (n-i-1 < delta) { // close to the end of the data vector
      d = n-i-1;
    } else { // somewhere in the middle of the data vector
      d = delta;
    }
    yy(i) = y.segment(i-d,2*d+1).sum()/(2*d+1); // perform filtering
  } // end iteration through the data

  return yy; // return new filtered data vector
}

// perform filtering by running average inplace
// static member function
void Curve::filterRunningAverage(VectorXd& y, const int span) {
  y = Curve::filteredRunningAverage(y, span); // get filtered data vector and
                                              // assign it to the old variable
}

// resample the vector using the number of points
// static member function
void Curve::resample(VectorXd& x, const int m,
                     const gsl_interp_type* interpType)
{
  x = resampled(x, m, interpType); // resample X inplace, by reference
}

// returns resampled vector using the number of points
// static member function
VectorXd Curve::resampled(const VectorXd& x, const int m,
                          const gsl_interp_type* interpType)
{
  const int n = x.size(); // vector length

  // chack if there are enough points for interpolation
  if (n >= gsl_interp_type_min_size(interpType)) {
    // initialize the accelerator object
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    // initialize the interpolating spline
    gsl_spline *spline = gsl_spline_alloc(interpType, n);

    const double *xArr, *tArr; // array for storing vector data

    xArr = x.data(); // save vector data to array
    // create linearly spaced vector - time
    VectorXd t = VectorXd::LinSpaced(n, 0, n-1);
    tArr = t.data(); // and save it to array

    // initialize the interpolating spline
    gsl_spline_init(spline, tArr, xArr, n);
    const double minT = t[0]; // starting time
    const double maxT = t[n-1]; // end time
    const double dt = (maxT-minT)/(m-1); // compute the dt step
    VectorXd xx = VectorXd::Zero(m); // initialize resampled vector object
    // save new resampled curve data
    for (int i = 0; i < m; i++) {
      // avoid extrapolation
      xx(i) = (i == m-1 ? gsl_spline_eval(spline, maxT, acc) :
                          gsl_spline_eval(spline, minT+i*dt, acc));
    }
    gsl_spline_free(spline); // free spline memory
    gsl_interp_accel_free(acc); // free accelerator memory

    return xx; // return resampled vector
  } else {
    // if interpolation requires more points
    throw "Not enough points for interpolation";
  }
}

// returns weighted average as in Hau & Sonnerup (1999)
// static member function
VectorXd Curve::weightedAverage(const VectorXd& x, double w) {
  // initialize weighted averaged vector
  VectorXd xx(x);

  // do the averaging
  xx.segment(1, xx.size()-2) = ((3-2*w)*xx.segment(1, xx.size()-2)+
    w*(xx.head(xx.size()-2)+xx.tail(xx.size()-2)))/3;
//  xx.segment(1, xx.size()-2) = w*xx.segment(1, xx.size()-2)+
//    0.5*(1-w)*(xx.head(xx.size()-2)+xx.tail(xx.size()-2));

  return xx; // return new averaged vector
}

// sort the vector inplace
// static member function
void Curve::sort(VectorXd& x) {
  x = sorted(x); // sort the vector inplace
}

// returns sorted vector
// static member function
VectorXd Curve::sorted(VectorXd& x) {

  double *xArr; // array for storing vector data
  xArr = x.data(); // write vector data to array

  gsl_sort(xArr, 1, x.size()); // sort the array
  Map<VectorXd> xx(xArr, x.size()); // initialize new sorted vector

  return xx; // return the sorted vector
}


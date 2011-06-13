
#ifndef CURVE_H
#define CURVE_H

#include <eigen3/Eigen/Dense>
#include <gsl/gsl_interp.h>

#include <vector>

// this structure stores curve data
struct CurveVectors {
  Eigen::VectorXd x, y; // x and y components of the curve
};

// this class is used for storing information for a single curve, which is
// basically (x,y) paired data, the curve is not a function, i.e. y may be not
// unique for a single x value
class Curve {
  protected:
    CurveVectors _vectors; // curve data
  public:
    Curve(); // constructor
    Curve(Eigen::VectorXd, Eigen::VectorXd); // construct out of vectors
    // curve data accessor
    const CurveVectors& cols() const {return _vectors;}
    const int size() const {return _vectors.x.size();} // length of the curve
    Curve& filterRunningAverage(int, char); // running average filter
    Curve& filterSavitzkyGolay(int, int); // Savitzky-Golay filter
    // resample the curve using min and max X limits and number of points
    Curve& resample(const double, const double, const int,
                    const gsl_interp_type* interpType = gsl_interp_linear);
    // resample the curve using only number of points required
    Curve& resample(const int,
                    const gsl_interp_type* interpType = gsl_interp_linear);
    // compute the residue between the branches of the curve
    double computeResidue();
    // returns data vector smoothed by running average filter
    static Eigen::VectorXd filteredRunningAverage(const Eigen::VectorXd&,
                                                  const int);
    // filters data vector by running average inplace
    static void filterRunningAverage(Eigen::VectorXd&, const int);
    // returns data vector smoothed by Savitzky-Golay filter
    static Eigen::VectorXd filteredSavitzkyGolay(const Eigen::VectorXd&,
                                                 const int, const int);
    // filters data vector by Savitzky-Golay filter
    static void filterSavitzkyGolay(Eigen::VectorXd&, const int, const int);
    // resample the vector using only number of points required
    static void resample(Eigen::VectorXd&, const int,
                    const gsl_interp_type* interpType = gsl_interp_linear);
    // returns resampled vector using only the number of points
    static Eigen::VectorXd resampled(const Eigen::VectorXd&, const int,
                    const gsl_interp_type* interpType = gsl_interp_linear);
    // returns weighted average as in Hau & Sonnerup (1999)
    static Eigen::VectorXd weightedAverage(const Eigen::VectorXd&, double);
};

#endif


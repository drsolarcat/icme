
#ifndef CURVE_H
#define CURVE_H

#include <eigen3/Eigen/Dense>

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
    CurveVectors c_vectors; // curve data
  public:
    Curve(); // constructor
    // curve data accessor
    const CurveVectors& cols() const {return c_vectors;}
    Curve& filterRunningAverage(int, char); // running average filter
    Curve& filterSavitzkyGolay(int, int); // Savitzky-Golay filter
    // compute the residue between the branches of the curve
    double computeResidue();
    // returns data vector smoothed by running average filter
    static Eigen::VectorXd filteredRunningAverage(Eigen::VectorXd, int);
    // filters data vector by running average inplace
    static void filterRunningAverage(Eigen::VectorXd&, int);
};

#endif


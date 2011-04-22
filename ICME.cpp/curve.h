
#ifndef CURVE_H
#define CURVE_H

#include <vector>

// this class is used for storing information for a single curve, which is
// basically (x,y) paired data, the curve is not a function, i.e. y may be not
// unique for a single x value
class Curve {
    vector<double> x, y;
  public:
    Curve(); // constructor
    int size(); // size of the curve, i.e. the number of elements in x and y
                // vectors
    Curve& smooth(); // smooth the curve, running average or Golay filter
    Curve copy(); // copy the curve object, to be used before smoothing
};

#endif


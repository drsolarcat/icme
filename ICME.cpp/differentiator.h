
#ifndef DIFFERENTIATOR
#define DIFFERENTIATOR

// library headers
#include <eigen3/Eigen/Dense>

class Differentitor {
  public:
    Eigen::VectorXd Holoborodko2(const int,
                                 const Eigen::VectorXd&, const double);
};

#endif


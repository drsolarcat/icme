
#ifndef DIFFERENTIATOR
#define DIFFERENTIATOR

// library headers
#include <eigen3/Eigen/Dense>

// a class for numerical differetitation
class Differentitor {
  public:
    // Holoborodko 2nd derivative filter
    Eigen::VectorXd Holoborodko2(const int,
                                 const Eigen::VectorXd&, const double);
};

#endif


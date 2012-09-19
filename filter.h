
#ifndef FILTER_H
#define FILTER_H

// library headers
#include <eigen3/Eigen/Dense>

// class for filtering data
class Filter {
  public:
    // 1D median filter, good for removing spikes from data
    static Eigen::VectorXd median1D(Eigen::VectorXd, int);
};

#endif


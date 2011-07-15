
#include "curve.h"
#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

int main() {

  const int n = 10;
  VectorXd y = VectorXd::Zero(n);

  VectorXd x(y);

  x.segment(0, 5) = Curve::weightedAverage(y, 0.5);

  cout << y << endl << endl << x << endl;

  return 0;
}


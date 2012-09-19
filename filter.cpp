
// project headers
#include "filter.h"
// standard headers
#include <vector>

using namespace Eigen;
using namespace std;

VectorXd Filter::median1D(VectorXd x, int window) {
  VectorXd y = x;
  if (window%2 == 0) window+=1;
  int span = ceil(window/2);
  for (int i = 0; i <x.size(); i++) {
    vector<double> v(window);
    for (int k = -span; k <= span; k++) {
      if (i+k < 0) {
        v[k+span] = x(0);
      } else if (i+k > x.size()-1) {
        v[k+span] = x(x.size()-1);
      } else {
        v[k+span] = x(i+k);
      }
    }
    sort(v.begin(), v.begin()+window);
    y(i) = v[span+1];
  }
  return y;
}


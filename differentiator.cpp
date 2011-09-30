
// project headers
#include "differentiator.h"
// standard headers
#include <iostream>

using namespace std;
using namespace Eigen;

// Holoborodko 2nd derivative
VectorXd Differentitor::Holoborodko2(const int span,
                                     const VectorXd& x, const double dx)
{
  const int n = x.size(); // length of the vector
  VectorXd res = VectorXd::Zero(n); // initialize the result

  int m = 0;
  VectorXd s;
  if (span%2 == 1) { // only for odd span
    for (int w = span; w >= 1; w = w-2) { // loop through spans
      if (w == 1) { // first and last points of the vector
        res(0)   = (2*x(0)-5*x(1)+4*x(2)-x(3))/pow(dx, 2); // first
        res(n-1) = (2*x(n-1)-5*x(n-2)+4*x(n-3)-x(n-4))/pow(dx, 2); // last
      } else {
        m = (w-1)/2; // center of span
        // initialize and fill helper coefficients
        s = VectorXd::Zero(m+2);
        s(m) = 1;
        for (int k = m-1; k >= 0; k--) {
          s(k) = ((2*w-10)*s(k+1)-(w+2*k+3)*s(k+2))/(w-2*k-1);
        }
        s.conservativeResize(m+1);
        if (w == span) { // central part
          for (int i = m; i <= n-1-m; i++) {
            res(i) = (s(0)*x(i)+((x.segment(i+1, m)+
              x.segment(i-m, m).reverse()).array()*s.tail(m).array()).sum())/
              pow(2,w-3)/pow(dx,2);
          }
        } else { // periferal parts
          res(m)     = (s(0)*x(m)+((x.segment(m+1, m)+
            x.segment(0, m).reverse()).array()*s.tail(m).array()).sum())/
            pow(2,w-3)/pow(dx,2);
          res(n-1-m) = (s(0)*x(n-1-m)+((x.segment(n-m, m)+
            x.segment(n-1-2*m, m).reverse()).array()*s.tail(m).array()).sum())/
            pow(2,w-3)/pow(dx,2);
        }
      }
    }
  }

  return res; // return the vector of 2nd derivatives
}


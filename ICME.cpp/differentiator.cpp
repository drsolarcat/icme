
// project headers
#include "differentiator.h"

using namespace Eigen;

VectorXd Differentitor::Holoborodko2(const int span,
                                     const VectorXd& x, const double dx)
{
  const int n = x.size();
  VectorXd res = VectorXd::Zero(n);

  int m = 0;
  VectorXd s;
  if (span%2 == 1) {
    for (int w = span; w >= 1; w = w-2) {
      if (w == 1) {
        res(0)   = (2*x(0)-5*x(1)+4*x(2)-x(3))/pow(dx, 2);
        res(n-1) = (2*x(n-1)-5*x(n-2)+4*x(n-3)-x(n-4))/pow(dx, 2);
      } else {
        m = (w-1)/2;
        s = VectorXd::Zero(m+2);
        for (int k = m-1; k >= 0; k--) {
          s(k) = ((2*w-10)*s(k+1)-(w+2*k+3)*s(k+2))/(w-2*k-1);
        }
        s.conservativeResize(m+1);
        if (w == span) {
          for (int i = m; i <= n-1-m; i++) {
            res(i) = s(0)*x(i)+((x.segment(i+1, m)+
              x.segment(i-m, m).reverse()).array()*s.tail(m).array()).sum();
          }
        } else {
          res(m)     = s(0)*x(m)+((x.segment(m+1, m)+
            x.segment(0, m).reverse()).array()*s.tail(m).array()).sum();
          res(n-1-m) = s(0)*x(n-1-m)+((x.segment(n-m, m)+
            x.segment(n-1-2*m, m).reverse()).array()*s.tail(m).array()).sum();
        }
      }
    }
  }

  return res;
}


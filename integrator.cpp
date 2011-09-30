
// project headers
#include "integrator.h"

using namespace Eigen;

Integrator::Integrator() {
  // trapezoid rule
  _NewtonCotesCoefficients[2].resize(2);
  _NewtonCotesCoefficients[2] << 1, 1;
  _NewtonCotesCoefficients[2] /= 2;
  // Simpson's rule
  _NewtonCotesCoefficients[3].resize(3);
  _NewtonCotesCoefficients[3] << 1, 4, 1;
  _NewtonCotesCoefficients[3] /= 3;
  // Simpson's 3/8 rule
  _NewtonCotesCoefficients[4].resize(4);
  _NewtonCotesCoefficients[4] << 1, 3, 3, 1;
  _NewtonCotesCoefficients[4] *= double(3)/double(8);
  // Boole's rule
  _NewtonCotesCoefficients[5].resize(5);
  _NewtonCotesCoefficients[5] << 7, 32, 12, 32, 7;
  _NewtonCotesCoefficients[5] *= double(2)/double(45);
  // Newton-Cotes N=6 rule
  _NewtonCotesCoefficients[6].resize(6);
  _NewtonCotesCoefficients[6] << 19, 75, 50, 50, 75, 19;
  _NewtonCotesCoefficients[6] *= double(5)/double(288);
  // Holoborodko's N=5 rule
  _HoloborodkoCoefficients[5].resize(5);
  _HoloborodkoCoefficients[5] << 11, 26, 31, 26, 11;
  _HoloborodkoCoefficients[5] *= double(4)/double(105);
  // Holoborodko's N=6 rule
  _HoloborodkoCoefficients[6].resize(6);
  _HoloborodkoCoefficients[6] << 31, 61, 76, 76, 61, 31;
  _HoloborodkoCoefficients[6] *= double(5)/double(336);
  // Holoborodko's N=7 rule
  _HoloborodkoCoefficients[7].resize(7);
  _HoloborodkoCoefficients[7] << 7, 12, 15, 16, 15, 12, 7;
  _HoloborodkoCoefficients[7] /= 14;
}

// Newton-Cotes integrator
double Integrator::NewtonCotes(int span,
                               const VectorXd& x,
                               const VectorXd& y)
{
  const double dx = x(1)-x(0); // x step size
  const int n = x.size(); // size of data
  const int segments = (int)floor((n-1)/(span-1)); // number of segments
  double res = 0; // result of integration

  // iterate though the segments
  for (int i = 0; i < segments; i++) {
    // apply Newton-Cotes formulas
    res += (_NewtonCotesCoefficients[span].array()*
            y.segment(i*(span-1), span).array()).sum();
  } // end iterating through the segments
  res *= dx; // multiply by segment width

  // how many points left out of segments
  const int m = n - segments*(span-1) - 1;

  // integrate the remained points with lower order rule
  if (m > 0) res += NewtonCotes(span-1, x.tail(m+1), y.tail(m+1));

  return res; // return the result of the integration
}

// Holoborodko integrator
double Integrator::Holoborodko(int span,
                               const VectorXd& x,
                               const VectorXd& y)
{
  const double dx = x(1)-x(0); // x step size
  const int n = x.size(); // size of data
  const int segments = (int)floor((n-1)/(span-1)); // number of segments
  double res = 0; // result of integration

  // iterate though the segments
  for (int i = 0; i < segments; i++) {
    // apply Holoborodko formulas
    res += (_HoloborodkoCoefficients[span].array()*
            y.segment(i*(span-1), span).array()).sum();
  } // end iterating through the segments
  res *= dx; // multiply by segment width

  // how many points left out of segments
  const int m = n - segments*(span-1) - 1;

  // integrate the remained points with lower order rule
  if (m > 0) res += NewtonCotes(span-1, x.tail(m+1), y.tail(m+1));

  return res; // return the result of the integration
}


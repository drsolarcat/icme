
// project headers
#include "fit_poly_exp.h"
#include "fit_poly.h"
#include "fit_exp.h"
// library headers
#include <eigen3/Eigen/Dense>
// standard headers
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;

// fit method
void PolyExpFit::fit()
{
  PolyFit lineFit(_n, _x, _y, 1);
  _slope = (lineFit.c()[1] > 0 ? 1 : -1);
  _xBdr = (_slope > 0 ? _x[0] : _x[_n-1]);
  _xCtr = (_slope > 0 ? _x[_n-1] : _x[0]);

  _polyFit = new PolyFit(_n, _x, _y, _order);
  _polyFit->fit();
  _expFitBdr = new PosZeroExpFit(_n, _x, _y);
  _expFitBdr->fit();
  const int n = 50;
  double mx = (_x[_n-1]-_x[0])/2;
  double xc[n], yc[n];
  double x1 = (_slope > 0 ? mx : _x[0]);
  double x2 = (_slope > 0 ? _x[_n-1] : mx);
  double dx = (x2-x1)/n;
  double x = x1;
  int i = 0;
  while (x <= x2) {
    xc[i] = x;
    yc[i] = _polyFit->f(x);
    x += dx;
    i++;
  }
  _expFitCtr = new PosExpFit(n, xc, yc);
  _expFitCtr->fit();
}

// evaluate function at point X
double PolyExpFit::f(double x)
{
  return _expFitBdr->f(x)/(1+exp(-_slope/_sBdr*(-x+_xBdr)))+
         _polyFit->f(x)/(1+exp(-_slope/_sBdr*(x-_xBdr)))/
                       (1+exp(-_slope/_sCtr*(-x+_xCtr)))+
         _expFitCtr->f(x)/(1+exp(-_slope/_sCtr*(x-_xCtr)));
}

// evaluate function derivative at X
double PolyExpFit::df(double x)
{
  return _expFitBdr->df(x)/(1+exp(-_slope/_sBdr*(-x+_xBdr)))-
         _expFitBdr->f(x)*exp(-_slope/_sBdr*(-x+_xBdr))*_slope/_sBdr/
         pow(1+exp(-_slope/_sBdr*(-x+_xBdr)), 2)+
         _polyFit->df(x)/(1+exp(-_slope/_sBdr*(x-_xBdr)))/
                         (1+exp(-_slope/_sCtr*(-x+_xCtr)))-
         _polyFit->f(x)*
         (-exp(-_slope/_sBdr*(x-_xBdr))*_slope/_sBdr*
         (1+exp(-_slope/_sCtr*(-x+_xCtr)))+
         exp(-_slope/_sCtr*(-x+_xCtr))*_slope/_sCtr*
         (1+exp(-_slope/_sBdr*(x-_xBdr))))/
         pow((1+exp(-_slope/_sBdr*(x-_xBdr)))*
             (1+exp(-_slope/_sCtr*(-x+_xCtr))), 2)+
         _expFitCtr->df(x)/(1+exp(-_slope/_sCtr*(x-_xCtr)))+
         _expFitCtr->f(x)*exp(-_slope/_sCtr*(x-_xCtr))*_slope/_sCtr/
         pow(1+exp(-_slope/_sCtr*(x-_xCtr)), 2);
}


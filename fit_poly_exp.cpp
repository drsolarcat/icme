
// project headers
#include "fit_poly_exp.h"
#include "fit_poly.h"
#include "fit_exp.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <log4cplus/logger.h>
#include <log4cplus/configurator.h>
// standard headers
#include <cmath>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace log4cplus;

// fit method
void PolyExpFit::fit()
{

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  // fit linearly to find the slope
  PolyFit lineFit(_n, _x, _y, 1);
  lineFit.fit();
  // slope
  _slope = (lineFit.c()[1] > 0 ? 1 : -1);
  LOG4CPLUS_TRACE(logger, "the curve slope = " << _slope);
  // boundaries of the curve
  _xBdr = (_slope > 0 ? _x[0] : _x[_n-1]);
  _xCtr = (_slope > 0 ? _x[_n-1] : _x[0]);
  LOG4CPLUS_TRACE(logger, "boundary x-value = " << _xBdr);
  LOG4CPLUS_TRACE(logger, "boundary s-param = " << _sBdr);
  LOG4CPLUS_TRACE(logger, "center x-value = " << _xCtr);
  LOG4CPLUS_TRACE(logger, "center s-param = " << _sCtr);

  // fit the the polynomial part of the curve
  _polyFit = new PolyFit(_n, _x, _y, _order);
  _polyFit->fit();
  LOG4CPLUS_TRACE(logger, "polynomial order = " << _order);
  if (_order == 2) {
    LOG4CPLUS_TRACE(logger, "polynomial params = [" <<
                            _polyFit->c()[0] << ", " <<
                            _polyFit->c()[1] << ", " <<
                            _polyFit->c()[2] << "]");
  } else if (_order == 3) {
    LOG4CPLUS_TRACE(logger, "polynomial params = [" <<
                            _polyFit->c()[0] << ", " <<
                            _polyFit->c()[1] << ", " <<
                            _polyFit->c()[2] << ", " <<
                            _polyFit->c()[3] << "]");
  }

  // fit the exponential tail at the boundary, i.e. lower end of the curve
  _expFitBdr = new PosZeroExpFit(_n, _x, _y);
  _expFitBdr->fit();
  LOG4CPLUS_TRACE(logger, "boundary exponent params = [" <<
                          _expFitBdr->c()[0] << ", " <<
                          _expFitBdr->c()[1] << "]");

  // fit the exponential tail at the center, i.e. upper end of the curve
  const int n = 50; // number of points to fit
  double mx = _x[0]+(_x[_n-1]-_x[0])/2; // center of the curve
  double xc[n], yc[n]; // new arrays for points to be fitted
  double x1 = (_slope > 0 ? mx : _x[0]); // left end of the fitting part
  double x2 = (_slope > 0 ? _x[_n-1] : mx); // right end of the fitting part
  LOG4CPLUS_DEBUG(logger, "x1 = " << x1 << ", x2 = " << x2);
  double dx = (x2-x1)/(n-1); // step size
  double x = x1; // initial position
  int i = 0;
  for (int i = 0; i < n; i++) { // fill in the temporary arrays to be fitted
    xc[i] = x1+i*dx;
    yc[i] = _polyFit->f(xc[i]);
  }
  // finally fit the exponential at the center
  _expFitCtr = new PosExpFit(n, xc, yc);
  _expFitCtr->fit();
  LOG4CPLUS_TRACE(logger, "center exponent params = [" <<
                          _expFitCtr->c()[0] << ", " <<
                          _expFitCtr->c()[1] << ", " <<
                          _expFitCtr->c()[2] << "]");
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


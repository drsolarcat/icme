
// project headers
#include "fit_exp.h"
#include "curve.h"
#include "integrator.h"
// library headers
#include <eigen3/Eigen/Dense>

using namespace Eigen;

void ExpFit::fit0()
{
  // map array data to vectors
  Map<VectorXd> xVec(_x, _n);
  Map<VectorXd> yVec(_y, _n);

  // compute minimum and maximum X values
  double xMin = xVec.minCoeff(),
         xMax = xVec.maxCoeff();
  // range of the curve in X
  double X = xMax-xMin;
  // resampling length
  int m = 200;

  // initialize the curve out of X and Y vectors
  Curve curve(xVec, yVec);
  // resample the curve
  curve.resample(m);

  // initialize integrator object
  Integrator integrator;
  // initialize coefficients
  double I1, I2, I3, I4;
  // and compute them
  I1 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(0)*curve.cols().y.array())/1;
  I2 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(1)*curve.cols().y.array())/1;
  I3 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(2)*curve.cols().y.array())/2;
  I4 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(3)*curve.cols().y.array())/6;

  // finally, compute the initial guess parameters
  double tau = (12*I4-6*I3*X+I2*pow(X, 2))/(-12*I3+6*I2*X-I1*pow(X, 2));
  double Q1 = exp(-xMin/tau);
  double Q = exp(-X/tau);
  _c0[0] = 2/X/((1+Q)*X+2*(Q-1)*tau)*(I2*(Q-1)+I1*(X+(Q-1)*tau));
  _c0[1] = (2*I2-I1*X)/tau/((1+Q)*X+2*(Q-1)*tau)/Q1;
  _c0[2] = -1/tau;
}

void PosExpFit::fit0()
{
  _c0[0] = 0;
  _c0[2] = log(_y[0]/_y[_n-1])/(_x[0]-_x[_n-1]);
  _c0[1] = log(_y[0])-_c0[2]*_x[0];
}

void PosZeroExpFit::fit0()
{
  _c0[1] = log(_y[0]/_y[_n-1])/(_x[0]-_x[_n-1]);
  _c0[0] = log(_y[0])-_c0[1]*_x[0];
}


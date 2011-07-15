
// project headers
#include "fit_cexp.h"
#include "curve.h"
#include "integrator.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
// standard headers
#include <iostream>

using namespace std;
using namespace Eigen;

// process the fitting
int fit_cexp(int n, double* x, double* y, double* c)
{
  const gsl_multifit_fdfsolver_type *T; // solver type
  gsl_multifit_fdfsolver *s; // solver
  int status;
  const int nc = 3; // number of fitting coefficients

  double c0[3]; // initialize initial guess for fitting parameters
  fit_cexp_0(n, x, y, c0); // compute the initial guess

  gsl_multifit_function_fdf f; // fitting function
  fit_cexp_data params; // fitting function parameters

  // initialize parametersof the fitting function
  params.x = x;
  params.y = y;
  params.n = n;

  // initialize the fitting function
  f.f = &fit_cexp_f; // function
  f.df = &fit_cexp_df; // Jacobian
  f.fdf = &fit_cexp_fdf; // function and Jacobian
  f.n = n; // length
  f.p = nc; // number of fitting coefficients
  f.params = &params; // fitting function parameters

  // vector for initial guess coefficients
  gsl_vector* cVec0 = gsl_vector_alloc(nc);
  gsl_vector_set(cVec0, 0, c0[0]);
  gsl_vector_set(cVec0, 1, c0[1]);
  gsl_vector_set(cVec0, 2, c0[2]);

  // set the solver
  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, nc);
  gsl_multifit_fdfsolver_set(s, &f, cVec0);

  // iterate the solver
  int i = 0;
  do {
    i++;
    status = gsl_multifit_fdfsolver_iterate(s);
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-40, 1e-40);
  } while (status == GSL_CONTINUE && i < 5000);

  // coefficients of the fitted curve
  c[0] = gsl_vector_get(s->x, 0);
  c[1] = gsl_vector_get(s->x, 1);
  c[2] = gsl_vector_get(s->x, 2);

  // free the solver
  gsl_multifit_fdfsolver_free (s);

  return GSL_SUCCESS; // done
}

// compute the initial guess for the fitting parameters
int fit_cexp_0(int n, double* x, double* y, double* c0)
{
  // map array data to vectors
  Map<VectorXd> xVec(x, n);
  Map<VectorXd> yVec(y, n);

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
  c0[0] = 2/X/((1+Q)*X+2*(Q-1)*tau)*(I2*(Q-1)+I1*(X+(Q-1)*tau));
  c0[1] = (2*I2-I1*X)/tau/((1+Q)*X+2*(Q-1)*tau)/Q1;
  c0[2] = -1/tau;

  return GSL_SUCCESS; // done
}

// minimization function
int fit_cexp_f(const gsl_vector* coeff, void* params, gsl_vector* f)
{
  // get function parameters
  double* x = ((fit_cexp_data*)params)->x;
  double* y = ((fit_cexp_data*)params)->y;
  int n = ((fit_cexp_data*)params)->n;

  const int nc = 3; // number of fitting coefficients
  double c[nc]; // coefficients array
  // fill the coefficients array
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double Yi; // fitting function value

  // fill the minimization function values
  for (int i = 0; i < n; i++) {
    Yi = c[0]+c[1]*exp(c[2]*x[i]);
    gsl_vector_set(f, i, Yi-y[i]);
  }

  return GSL_SUCCESS; // done
}

// compute Jacobian of the fitting function
int fit_cexp_df(const gsl_vector* coeff, void* params, gsl_matrix* J)
{
  // get function parameters from the parameters structure
  double* x = ((fit_cexp_data*)params)->x;
  double* y = ((fit_cexp_data*)params)->y;
  int n = ((fit_cexp_data*)params)->n;

  const int nc = 3; // number of  fitting coefficients
  double c[nc]; // initialize array of the fitting coefficients
  // fill the coefficients array
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double df_dc0, df_dc1, df_dc2; // derivatives over fitting coefficients
  // fill the Jacobian with derivatives
  for (int i = 0; i < n; i++) {
    df_dc0 = 1; // df/dc0
    df_dc1 = exp(c[2]*x[i]); // df/dc1
    df_dc2 = c[1]*x[i]*exp(c[2]*x[i]); // df/dc2

    gsl_matrix_set(J, i, 0, df_dc0);
    gsl_matrix_set(J, i, 1, df_dc1);
    gsl_matrix_set(J, i, 2, df_dc2);
  }

  return GSL_SUCCESS; // done
}

// compute Jacobian and function simultaiously
int fit_cexp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                 gsl_matrix* J)
{
  fit_cexp_f(coeff, params, f); // compute the minimization function
  fit_cexp_df(coeff, params, J); // compute the Jacobian

  return GSL_SUCCESS; // done
}

// evaluate the fitting function
double fit_cexp_eval_f(double x, double* c)
{
  return c[0]+c[1]*exp(c[2]*x); // return the f(x) value
}

// evaluate the derivative of the fitting function
double fit_cexp_eval_df(double x, double* c)
{
  return c[1]*c[2]*exp(c[2]*x); // return the df/dx(x) value
}


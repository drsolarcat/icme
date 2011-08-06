
// project headers
#include "fit_poly.h"
// library headers
#include <gsl/gsl_multifit.h>

// fit method
void PolyFit::fit()
{
  gsl_multifit_linear_workspace *ws; // workspace
  gsl_matrix *cov, *X; // covariance and X data matrices
  gsl_vector *yVec, *cVec; // Y data and coefficients vectors
  double chisq; // chi squared

  // allocate vectors and matrices
  X = gsl_matrix_alloc(_n, _nc);
  yVec = gsl_vector_alloc(_n);
  cVec = gsl_vector_alloc(_nc);
  cov = gsl_matrix_alloc(_nc, _nc);

  // fill X data matrix
  for(int i = 0; i < _n; i++) {
    gsl_matrix_set(X, i, 0, 1.0);
    for(int j = 0; j <= _order; j++) {
      gsl_matrix_set(X, i, j, pow(_x[i], j));
    }
    gsl_vector_set(yVec, i, _y[i]);
  }

  // allocate linear multifit solver
  ws = gsl_multifit_linear_alloc(_n, _nc);
  // do the fitting
  gsl_multifit_linear(X, yVec, cVec, cov, &chisq, ws);

  // fill in fitting coefficients array
  for(int i = 0; i <= _order; i++) {
    _c[i] = gsl_vector_get(cVec, i);
  }

  // free all preallocated objects
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(yVec);
  gsl_vector_free(cVec);
}

double PolyFit::f(double x, int order, double* c)
{
  double y = 0; // initialize the resulting value

  // sum the polynomial parts
  for (int i = 0; i <= order; i++) {
    y += c[i]*pow(x, i);
  }

  return y; // return the result
}

double PolyFit::f(double x)
{
  return f(x, _order, _c); // return the result
}

// evaluate function derivative at X
double PolyFit::df(double x)
{
  double cd[_nc-1];
  for (int i = 0; i < _nc-1; i++) {
    cd[i] = _c[i+1]*(i+1);
  }

  return f(x, _order-1, cd);
}



// project headers
#include "fit_poly.h"
// library headers
#include <gsl/gsl_multifit.h>

// fit polynomial of required order to data
int fit_poly(int n, double* x, double* y, int order, double* c)
{
  gsl_multifit_linear_workspace *ws; // workspace
  gsl_matrix *cov, *X; // covariance and X data matrices
  gsl_vector *yVec, *cVec; // Y data and coefficients vectors
  double chisq; // chi squared

  // allocate vectors and matrices
  X = gsl_matrix_alloc(n, order+1);
  yVec = gsl_vector_alloc(n);
  cVec = gsl_vector_alloc(order+1);
  cov = gsl_matrix_alloc(order+1, order+1);

  // fill X data matrix
  for(int i = 0; i < n; i++) {
    gsl_matrix_set(X, i, 0, 1.0);
    for(int j = 0; j <= order; j++) {
      gsl_matrix_set(X, i, j, pow(x[i], j));
    }
    gsl_vector_set(yVec, i, y[i]);
  }

  // allocate linear multifit solver
  ws = gsl_multifit_linear_alloc(n, order+1);
  // do the fitting
  gsl_multifit_linear(X, yVec, cVec, cov, &chisq, ws);

  // fill in fitting coefficients array
  for(int i = 0; i <= order; i++) {
    c[i] = gsl_vector_get(cVec, i);
  }

  // free all preallocated objects
  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(yVec);
  gsl_vector_free(cVec);

  return GSL_SUCCESS; // done
}

// evaluate polynomial
double fit_poly_eval_f(double x, int order, double* c)
{
  double y = 0; // initialize the resulting value

  // sum the polynomial parts
  for (int i = 0; i <= order; i++) {
    y += c[i]*pow(x, i);
  }

  return y; // return the result
}

// evaluate polynomial 1st derivative
double fit_poly_eval_df(double x, int order, double* c)
{
  double cd[order];
  for (int i = 0; i < order; i++) {
    cd[i] = c[i+1]*(i+1);
  }

  return fit_poly_eval_f(x, order-1, cd);
}


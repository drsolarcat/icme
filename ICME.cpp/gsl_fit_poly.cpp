
// project headers
#include "gsl_fit_poly.h"
// library headers
#include <gsl/gsl_multifit.h>

void gsl_fit_poly(const int n, const int order,
                  const double *xdata, const double *ydata,
                  double *coeff)
{
  gsl_multifit_linear_workspace *ws;
  gsl_matrix *cov, *X;
  gsl_vector *y, *c;
  double chisq;

  int i, j;

  X = gsl_matrix_alloc(n, order+1);
  y = gsl_vector_alloc(n);
  c = gsl_vector_alloc(order+1);
  cov = gsl_matrix_alloc(order+1, order+1);

  for(int i = 0; i < n; i++) {
    gsl_matrix_set(X, i, 0, 1.0);
    for(int j = 0; j <= order; j++) {
      gsl_matrix_set(X, i, j, pow(xdata[i], j));
    }
    gsl_vector_set(y, i, ydata[i]);
  }

  ws = gsl_multifit_linear_alloc(n, order+1);
  gsl_multifit_linear(X, y, c, cov, &chisq, ws);

  for(int i = 0; i <= order; i++)
  {
    coeff[i] = gsl_vector_get(c, i);
  }

  gsl_multifit_linear_free(ws);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(y);
  gsl_vector_free(c);
}

double gsl_fit_poly_eval(const double x, const double *coeff, const int order) {
  double y = 0;

  for (int i = 0; i <= order; i++) {
    y += coeff[i]*pow(x, i);
  }

  return y;
}


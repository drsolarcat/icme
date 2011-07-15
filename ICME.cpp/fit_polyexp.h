
#ifndef FIT_POLYEXP
#define FIT_POLYEXP

// project headers
// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// structure for storing the curve data and parameters
struct fit_polyexp_data {
  double *x, *y; // data
  int order;
  int n;
};

int fit_polyexp(int n, double* x, double* y, int order, double* c);
double fit_polyexp_eval_f(double x, int order, double* c);
double fit_polyexp_eval_df(double x, int order, double* c);

#endif


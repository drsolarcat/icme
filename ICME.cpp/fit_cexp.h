
#ifndef FIT_CEXP
#define FIT_CEXP

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// structure for storing the curve data
struct fit_cexp_data {
  double *x, *y; // X and Y data
  int n; // length of the curve
};

// fit the curve
int fit_cexp(int n, double* x, double* y, double* c);
// compute initial guess
int fit_cexp_0(int n, double* x, double* y, double* c0);
// fitting function
int fit_cexp_f(const gsl_vector* coeff, void* params, gsl_vector* f);
// compute Jacobian
int fit_cexp_df(const gsl_vector* coeff, void* params, gsl_matrix* J);
// compute fitting function and Jacobian
int fit_cexp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                       gsl_matrix* J);
// evaluate fitting function
double fit_cexp_eval_f(double x, double* c);
// evaluate the first derivative of the fitting function
double fit_cexp_eval_df(double x, double* c);

#endif


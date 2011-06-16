
#ifndef GSL_FIT_EXP
#define GSL_FIT_EXP

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// structure for storing the curve data
struct gsl_fit_exp_data {
  double *x, *y; // X and Y data
  int n; // length of the curve
};

// fit the curve
int gsl_fit_exp(int n, double* x, double* y, double* c);
// compute initial guess
int gsl_fit_exp_0(int n, double* x, double* y, double* c0);
// fitting function
int gsl_fit_exp_f(const gsl_vector* coeff, void* params, gsl_vector* f);
// compute Jacobian
int gsl_fit_exp_df(const gsl_vector* coeff, void* params, gsl_matrix* J);
// compute fitting function and Jacobian
int gsl_fit_exp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                    gsl_matrix* J);
// evaluate fitting function
double gsl_fit_exp_eval_f(double x, double* c);
// evaluate the first derivative of the fitting function
double gsl_fit_exp_eval_df(double x, double* c);

#endif


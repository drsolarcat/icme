
#ifndef GSL_FIT_EXP
#define GSL_FIT_EXP

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

struct gsl_fit_exp_data {
  double *x, *y;
  int n;
};

int gsl_fit_exp(int n, double* x, double* y, double* c);
int gsl_fit_exp_0(int n, double* x, double* y, double* c0);
int gsl_fit_exp_f(const gsl_vector* coeff, void* params, gsl_vector* f);
int gsl_fit_exp_df(const gsl_vector* coeff, void* params, gsl_matrix* J);
int gsl_fit_exp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                    gsl_matrix* J);
double gsl_fit_exp_eval_f(double x, double* c);
double gsl_fit_exp_eval_df(double x, double* c);

#endif


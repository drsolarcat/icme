
#ifndef GSL_FIT_EPE
#define GSL_FIT_EPE

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

struct gsl_fit_epe_data {
  double *x, *y;
  double *p, *e, *E;
  int order;
  double xc, xb;
  int n;
};

int gsl_fit_epe(int n, double* x, double* y, double xc, double xb,
                double* p, int order, double* e, double* E,
                double* c0, double* c);
int gsl_fit_epe_f(const gsl_vector* coeff, void* params, gsl_vector* f);
int gsl_fit_epe_df(const gsl_vector* coeff, void* params, gsl_matrix* J);
int gsl_fit_epe_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                    gsl_matrix* J);
double gsl_fit_epe_eval_f(double x, double xc, double xb, double* c, double* p,
                          int order, double* e, double* E);
double gsl_fit_epe_eval_df(double x, double xc, double xb, double* c,
                           double* p, double *dp, int order,
                           double* e, double* E);

#endif


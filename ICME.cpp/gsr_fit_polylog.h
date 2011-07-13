
#ifndef GSR_FIT_POLYLOG
#define GSR_FIT_POLYLOG

// project headers
// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// structure for storing the curve data and parameters
struct gsr_fit_polylog_data {
  double *x, *y; // data
  int order;
  int n;
};

int gsr_fit_polylog(const int n, const double* x, const double* y,
                    const int order, double* c);
int gsr_fit_polylog_0(const int n, const double* x, const double* y,
                      const int order, double* c0);
int gsl_fit_polylog_f(const gsl_vector* coeff, void* params, gsl_vector* f);
int gsl_fit_polylog_df(const gsl_vector* coeff, void* params, gsl_matrix* J);
int gsl_fit_polylog_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                        gsl_matrix* J);
double gsl_fit_polylog_eval_f(double x, const double* c, const int order);
double gsl_fit_polylog_eval_df(double x, const double* c, const int order);

#endif



#ifndef FIT_EXP
#define FIT_EXP

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// structure for storing the curve data
struct fit_exp_data {
  double *x, *y; // X and Y data
  int n; // length of the curve
};

// fit the curve
int fit_exp(int n, double* x, double* y, double* c);
// compute initial guess
int fit_exp_0(int n, double* x, double* y, double* c0);
// fitting function
int fit_exp_f(const gsl_vector* coeff, void* params, gsl_vector* f);
// compute Jacobian
int fit_exp_df(const gsl_vector* coeff, void* params, gsl_matrix* J);
// compute fitting function and Jacobian
int fit_exp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                gsl_matrix* J);
// evaluate fitting function
double fit_exp_eval_f(double x, double* c);
// evaluate the first derivative of the fitting function
double fit_exp_eval_df(double x, double* c);

class ExpFit {
  double *_c;
  public:
    ExpFit(int n, double *x, double *y) {
      _c = new double[n];
      fit_exp(n, x, y, _c);
    }
    ~ExpFit() {delete [] _c;}
    double f(double x) {return fit_exp_eval_f(x, _c);}
    double df(double x) {return fit_exp_eval_df(x, _c);}
};

#endif


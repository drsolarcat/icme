
#ifndef FIT_EXP
#define FIT_EXP

// library headers
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>

class ExpFit : public AbstractFit {
  protected:
    int _nc = 3;
  public:
    // evaluate function at point X using free parameters in array C
    virtual double f(double x, double* c);
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    virtual void df(double x, double* c, double* dy);
    // evaluate function derivative at X
    virtual double df(double x);
};

class PosExpFit : public AbstractFit {
  protected:
    int _nc = 3;
  public:
    // evaluate function at point X using free parameters in array C
    virtual double f(double x, double* c);
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    virtual void df(double x, double* c, double* dy);
    // evaluate function derivative at X
    virtual double df(double x);
};

class PosZeroExpFit : public AbstractFit {
  protected:
    int _nc = 2;
  public:
    // evaluate function at point X using free parameters in array C
    virtual double f(double x, double* c);
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    virtual void df(double x, double* c, double* dy);
    // evaluate function derivative at X
    virtual double df(double x);
};

/*
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
*/

#endif


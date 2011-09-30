
#ifndef FIT_ABSTRACT
#define FIT_ABSTRACT

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
// standard headers
#include <iostream>

// abstract class for fitting data to an analytical function
class AbstractFit {
  protected:
    int _n; // number of data points to fit
    int _nc; // number of free parameters to fit
    double *_x, *_y; // X and Y data
    int _niter; // number of iterations
    double _epsabs, _epsrel; // errors for stopping the fitting
    double *_c; // fitted free parameters
    double *_c0; // initial guess for free parameters
  public:
    // constructor
    AbstractFit() {}
    AbstractFit(int n, double* x, double* y, int nc, int niter,
                double epsabs, double epsrel) :
                _n(n), _x(x), _y(y), _nc(nc), _niter(niter),
                _epsabs(epsabs), _epsrel(epsrel)
    {
      // initialize coefficients arrays
      _c = new double[_nc];
      _c0 = new double[_nc];
    }
    // destructor
//    ~AbstractFit() {delete [] _c; delete [] _c0;}
    // fit method
    virtual void fit();
    // initial guess
    virtual void fit0() {for (int i = 0; i < _nc; i++) _c0[i] = 1;};
    // evaluate function at point X using free parameters in array C
    virtual double f(double x, double* c) = 0;
    // evaluate function at X
    double f(double x) {return f(x, _c);}
    // evaluate function at all points in array X
    void f(int n, double* x, double* y);
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    virtual void df(double x, double* c, double* dy) = 0;
    // evaluate function derivative at X with given coefficients C
    virtual double df(double x, double* c) = 0;
    // evaluate function derivative at X
    double df(double x) {return df(x, _c);}
    // evaluate function derivative at all points in array X
    void df(int n, double* x, double* dy);
    // accessors
    int n() const {return _n;}
    int nc() const {return _nc;}
    double* x() const {return _x;}
    double* y() const {return _y;}
    double* c() const {return _c;}
    // calculate vector of residues between fitted function and data
    int f(const gsl_vector* cVec, void* params, gsl_vector* fVec);
    // calculate Jacobian
    int df(const gsl_vector* cVec, void* params, gsl_matrix* J);
    // calculate both residues and Jacobian
    int fdf(const gsl_vector* cVec, void* params,
            gsl_vector* fVec, gsl_matrix* J);
};

#endif



#ifndef FIT_EXP
#define FIT_EXP

// project headers
#include "fit_abstract.h"
// standard headers
#include <cmath>

using namespace std;

class ExpFit : public AbstractFit {
  public:
    using AbstractFit::f;
    using AbstractFit::df;
    // constructor
    ExpFit(int n, double* x, double* y, int niter = 5000,
           double epsabs = 1e-40, double epsrel = 1e-40) :
           AbstractFit(n, x, y, 3, niter, epsabs, epsrel) {}
    // initial guess
    void fit0();
//    double f(double x) {return AbstractFit::f(x);}
    // evaluate function at point X using free parameters in array C
    double f(double x, double* c) {return c[0]+c[1]*exp(c[2]*x);}
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    void df(double x, double* c, double* dy) {
      dy[0] = 1;
      dy[1] = exp(c[2]*x);
      dy[2] = c[1]*x*exp(c[2]*x);
    }
    // evaluate function derivative at X
    double df(double x, double* c) {return c[1]*c[2]*exp(c[2]*x);}
};

class PosExpFit : public AbstractFit {
  public:
    using AbstractFit::f;
    using AbstractFit::df;
    // constructor
    PosExpFit(int n, double* x, double* y, int niter = 5000,
              double epsabs = 1e-40, double epsrel = 1e-40) :
              AbstractFit(n, x, y, 3, niter, epsabs, epsrel) {}
    // initial guess
    void fit0();
    // evaluate function at point X using free parameters in array C
    double f(double x, double* c) {return c[0]+exp(c[1]+c[2]*x);}
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    void df(double x, double* c, double* dy) {
      dy[0] = 1;
      dy[1] = exp(c[1]+c[2]*x);
      dy[2] = x*exp(c[1]+c[2]*x);
    }
    // evaluate function derivative at X
    double df(double x, double* c) {return c[2]*exp(c[1]+c[2]*x);}
};

class PosZeroExpFit : public AbstractFit {
  public:
    using AbstractFit::f;
    using AbstractFit::df;
    // constructor
    PosZeroExpFit(int n, double* x, double* y, int niter = 5000,
                  double epsabs = 1e-40, double epsrel = 1e-40) :
                  AbstractFit(n, x, y, 2, niter, epsabs, epsrel) {}
    // initial guess
    void fit0();
    // evaluate function at point X using free parameters in array C
    double f(double x, double* c) {return exp(c[0]+c[1]*x);}
    // evaluate function derivaives for all free parameters at point X using
    // free parameters in array C
    void df(double x, double* c, double* dy) {
      dy[0] = exp(c[0]+c[1]*x);
      dy[1] = x*exp(c[0]+c[1]*x);
    }
    // evaluate function derivative at X
    double df(double x, double* c) {return c[1]*exp(c[0]+c[1]*x);}
};

#endif


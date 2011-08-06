
#ifndef FIT_POLY_EXP
#define FIT_POLY_EXP

// project headers
#include "fit_poly.h"
#include "fit_exp.h"

class PolyExpFit {
  protected:
    int _n;
    double *_x, *_y;
    int _order;
    double _sCtr, _sBdr;
    double _xCtr, _xBdr;
    int _slope;
    PolyFit *_polyFit;
    PosZeroExpFit *_expFitBdr;
    PosExpFit *_expFitCtr;
  public:
    // constructor
    PolyExpFit(int n, double* x, double* y, int order,
               double sCtr, double sBdr) :
               _n(n), _x(x), _y(y), _order(order), _sCtr(sCtr), _sBdr(sBdr) {}
    // destructor
//    ~PolyExpFit() {delete _polyFit; delete _expFitBdr; delete _expFitCtr;}
    // fit method
    void fit();
    // evaluate function at point X
    double f(double x);
    // evaluate function derivative at X
    double df(double x);
};

#endif


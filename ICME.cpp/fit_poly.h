
#ifndef FIT_POLY
#define FIT_POLY

class PolyFit {
  protected:
    int _n;
    double *_x, *_y;
    int _order;
    int _nc;
    double *_c;
  public:
    // constructor
    PolyFit(int n, double* x, double* y, int order) :
      _n(n), _x(x), _y(y), _order(order), _nc(order+1)
    {
      _c = new double[_nc];
    }
    // fit method
    void fit();
    // evaluate function at point X using free parameters in array C
    double f(double x, int order, double* c);
    // evaluate function at point X
    double f(double x);
    // evaluate function derivative at X
    double df(double x);
    // coefficients accessor
    double* c() const {return _c;}
};

#endif


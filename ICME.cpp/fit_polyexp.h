
#ifndef FIT_POLYEXP
#define FIT_POLYEXP

// project headers
// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

// structure for storing the curve data and parameters
struct fit_polyexp_data {
  double *x, *y; // data
  int order;
  int n;
};

int fit_polyexp(int n, double* x, double* y, int order, double* c);
double fit_polyexp_eval_f(double x, int order, double* c);
double fit_polyexp_eval_df(double x, int order, double* c);

class PolyexpFit {
  double *_c;
  int _order;
  public:
    PolyexpFit(int n, double *x, double *y, int order) {
      _order = order;
      _c = new double[n];
      fit_polyexp(n, x, y, order, _c);
    }
    ~PolyexpFit() {delete [] _c;}
    double f(double x) {return fit_polyexp_eval_f(x, _order, _c);};
    double df(double x) {return fit_polyexp_eval_df(x, _order, _c);};
};

#endif


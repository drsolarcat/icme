
#ifndef GSL_FIT_POLY
#define GSL_FIT_POLY

void gsl_fit_poly(const int n, const int order,
                  const double *xdata, const double *ydata,
                  double *coeff);
double gsl_fit_poly_eval(const double x, const double *coeff, const int order);

#endif


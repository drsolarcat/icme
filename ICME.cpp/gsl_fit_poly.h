
#ifndef GSL_FIT_POLY
#define GSL_FIT_POLY

// fit polynomial
int gsl_fit_poly(const int n, const int order,
                 const double *xdata, const double *ydata,
                 double *coeff);
// evaluate polynomial
double gsl_fit_poly_eval(const double x, const double *coeff, const int order);

#endif


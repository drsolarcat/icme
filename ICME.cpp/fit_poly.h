
#ifndef FIT_POLY
#define FIT_POLY

// fit polynomial
int fit_poly(const int n, double* x, double* y, int order, double* c);
// evaluate polynomial
double fit_poly_eval_f(double x, int order, double* c);
// evaluate polynomial 1st derivative
double fit_poly_eval_df(double x, int order, double* c);
#endif


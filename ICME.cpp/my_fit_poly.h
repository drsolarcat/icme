
#ifndef MY_FIT_POLY
#define MY_FIT_POLY

// fit polynomial
int my_fit_poly(const int n, const double *x, const double *y,
                const int order, double *c);
// evaluate polynomial
double my_fit_poly_eval(const double x, const int order,
                        const double *c);

#endif


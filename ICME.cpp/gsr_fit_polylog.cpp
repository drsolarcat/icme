
// project headers
#include "gsr_fit_polylog.h"
#include "gsr_fit_poly.h"
#include "gsr_fit_exp.h"

int gsr_fit_polyexp(const int n, const double* x, const double* y,
                    const int order, double* c)
{
  double *cExp;
  gsr_fit_exp(n, x, y, order, cExp);
  double *cPoly;
  gsr_fit_poly(n, x, y, order, cPoly);
  for (int i = 0; i < order; i++) {
    c[i] = cPoly[i];
  }
  c[order] = cExp[0];
  c[order+1] = cExp[1];
  c[order+2] = 1;
  c[order+3] = x[0];
}

double gsl_fit_polyexp_eval_f(double x, const int order, const double* c)
{
  double *cPoly, *cExp, *cLog;
  for (int i = 0; i < order; i++) {
    cPoly[i] = c[i];
  }
  cExp[0] = c[order];
  cExp[1] = c[order+1];
  cLog[0] = c[order+2];
  cLog[1] = c[order+3];
  return gsr_fit_poly_eval(x, order, cPoly)/(1+exp(-cLog[0]*(x-cLog[1])))+
         gsr_fit_exp_eval(x, cExp)/(1+exp(-cLog[0]*(-x+cLog[1])));
}

double gsl_fit_polyexp_eval_df(double x, const int order, const double* c)
{

}


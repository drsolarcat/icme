
// project headers
#include "fit_polyexp.h"
#include "fit_poly.h"
#include "fit_exp.h"
// standard headers
#include <cmath>

int fit_polyexp(int n, double* x, double* y, int order, double* c)
{
  double cExp[2];
  fit_exp(n, x, y, cExp);
  double cPoly[order+1];
  fit_poly(n, x, y, order, cPoly);
  for (int i = 0; i < order+1; i++) {
    c[i] = cPoly[i];
  }
  c[order+1] = cExp[0];
  c[order+2] = cExp[1];
  c[order+3] = 0.08;
  c[order+4] = x[0];

  return GSL_SUCCESS;
}

double fit_polyexp_eval_f(double x, int order, double* c)
{
  double cPoly[order+1], cExp[2], cLog[2];

  for (int i = 0; i < order+1; i++) {
    cPoly[i] = c[i];
  }
  cExp[0] = c[order+1];
  cExp[1] = c[order+2];
  cLog[0] = c[order+3];
  cLog[1] = c[order+4];

  return fit_poly_eval_f(x, order, cPoly)/(1+exp(-cLog[0]*(x-cLog[1])))+
         fit_exp_eval_f(x, cExp)/(1+exp(-cLog[0]*(-x+cLog[1])));
}

double fit_polyexp_eval_df(double x, int order, double* c)
{
  double cPoly[order+1], cExp[2], cLog[2];

  for (int i = 0; i < order+1; i++) {
    cPoly[i] = c[i];
  }
  cExp[0] = c[order+1];
  cExp[1] = c[order+2];
  cLog[0] = c[order+3];
  cLog[1] = c[order+4];

  return (fit_poly_eval_df(x, order, cPoly)*(1+exp(-cLog[0]*(x-cLog[1])))+
         cLog[0]*exp(-cLog[0]*(x-cLog[1]))*fit_poly_eval_f(x, order, cPoly))/
         pow(1+exp(-cLog[0]*(x-cLog[1])),2)+
         (fit_exp_eval_df(x, cExp)*(1+exp(-cLog[0]*(-x+cLog[1])))-
         fit_exp_eval_f(x, cExp)*cLog[0]*exp(-cLog[0]*(-x+cLog[1])))/
         pow(1+exp(-cLog[0]*(-x+cLog[1])),2);

}


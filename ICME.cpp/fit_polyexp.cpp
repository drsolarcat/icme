
// project headers
#include "fit_polyexp.h"
#include "fit_poly.h"
#include "fit_exp.h"
// library headers
#include <gsl/gsl_sort.h>
#include <gsl/gsl_permute.h>
// standard headers
#include <cmath>

int fit_polyexp(int n, double* x, double* y, int order, double* c)
{
  // sort data
  size_t p[n];
  gsl_sort_index(p, x, 1, n);
  gsl_permute(p, x, 1, n);
  gsl_permute(p, y, 1, n);
  // fit polynomial part
  double cPoly[order+1];
  fit_poly(n, x, y, order, cPoly);
  // fit exponential part
  double cExp[2], cExpCtr[3], cExpBdr[2];
  fit_exp(n, x, y, cExp);

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


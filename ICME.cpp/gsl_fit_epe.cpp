
// project headers
#include "gsl_fit_epe.h"
#include "gsl_fit_poly.h"
// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
// standard headers
#include <iostream>
#include <limits>

using namespace std;

int gsl_fit_epe(int n, double* x, double* y, double xc, double xb,
                double* p, int order, double* e, double* E,
                double* c0, double* c)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  const int nc = 4;

  gsl_multifit_function_fdf f;
  gsl_fit_epe_data params;

  params.x = x;
  params.y = y;
  params.p = p;
  params.e = e;
  params.E = E;
  params.order = order;
  params.xc = xc;
  params.xb = xb;
  params.n = n;

  f.f = &gsl_fit_epe_f;
  f.df = &gsl_fit_epe_df;
  f.fdf = &gsl_fit_epe_fdf;
  f.n = n;
  f.p = nc;
  f.params = &params;

  gsl_vector *coeff0 = gsl_vector_alloc(nc);
  gsl_vector_set(coeff0, 0, c0[0]);
  gsl_vector_set(coeff0, 1, c0[1]);
  gsl_vector_set(coeff0, 2, c0[2]);
  gsl_vector_set(coeff0, 3, c0[3]);

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, nc);
  gsl_multifit_fdfsolver_set(s, &f, coeff0);

  int i = 0;
  do {
    i++;
    status = gsl_multifit_fdfsolver_iterate(s);
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE && i < 5000);

  c[0] = gsl_vector_get(s->x, 0);
  c[1] = gsl_vector_get(s->x, 1);
  c[2] = gsl_vector_get(s->x, 2);
  c[3] = gsl_vector_get(s->x, 3);

  gsl_multifit_fdfsolver_free (s);
}

int gsl_fit_epe_f(const gsl_vector* coeff, void* params, gsl_vector* f)
{
  double *x = ((gsl_fit_epe_data*)params)->x;
  double *y = ((gsl_fit_epe_data*)params)->y;
  double *p = ((gsl_fit_epe_data*)params)->p;
  double *e = ((gsl_fit_epe_data*)params)->e;
  double *E = ((gsl_fit_epe_data*)params)->E;
  int order = ((gsl_fit_epe_data*)params)->order;
  double xc = ((gsl_fit_epe_data*)params)->xc;
  double xb = ((gsl_fit_epe_data*)params)->xb;
  int n = ((gsl_fit_epe_data*)params)->n;

  const int nc = 4;
  double c[nc];
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double Yi;

  for (int i = 0; i < n; i++) {
    Yi = gsl_fit_poly_eval(x[i], p, order)/
         (1+c[0]*exp(c[1]*(x[i]-xb)))/
         (1+c[2]*exp(c[3]*(x[i]-xc)))+
         e[0]*exp(e[1]*x[i])/(1+c[0]*exp(c[1]*(-x[i]+xb)))+
         (E[0]*exp(E[1]*x[i])+E[2])/(1+c[2]*exp(c[3]*(-x[i]+xc)));
    gsl_vector_set(f, i, Yi-y[i]);
  }

  return GSL_SUCCESS;
}

int gsl_fit_epe_df(const gsl_vector* coeff, void* params, gsl_matrix* J)
{
  double *x = ((gsl_fit_epe_data*)params)->x;
  double *y = ((gsl_fit_epe_data*)params)->y;
  double *p = ((gsl_fit_epe_data*)params)->p;
  double *e = ((gsl_fit_epe_data*)params)->e;
  double *E = ((gsl_fit_epe_data*)params)->E;
  int order = ((gsl_fit_epe_data*)params)->order;
  double xc = ((gsl_fit_epe_data*)params)->xc;
  double xb = ((gsl_fit_epe_data*)params)->xb;
  int n = ((gsl_fit_epe_data*)params)->n;

  const int nc = 4;
  double c[nc];
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double df_dc0, df_dc1, df_dc2, df_dc3;
  for (int i = 0; i < n; i++) {
    df_dc0 = -gsl_fit_poly_eval(x[i], p, order)*exp(c[1]*(x[i]-xb))/
             pow(1+c[0]*exp(c[1]*(x[i]-xb)), 2)/
             (1+c[2]*exp(c[3]*(x[i]-xc)))
             -e[0]*exp(e[1]*x[i])*exp(c[1]*(-x[i]+xb))/
             pow(1+c[0]*exp(c[1]*(-x[i]+xb)), 2);
    df_dc1 = -gsl_fit_poly_eval(x[i], p, order)*exp(c[3]*(x[i]-xc))/
             (1+c[0]*exp(c[1]*(x[i]-xb)))/
             pow(1+c[2]*exp(c[3]*(x[i]-xc)), 2)
             -(E[0]*exp(E[1]*x[i])+E[2])*exp(c[3]*(-x[i]+xc))/
             pow(1+c[2]*exp(c[3]*(-x[i]+xc)), 2);
    df_dc2 = -gsl_fit_poly_eval(x[i], p, order)*
             c[0]*exp(c[1]*(x[i]-xb))*(x[i]-xb)/
             pow(1+c[0]*exp(c[1]*(x[i]-xb)), 2)/
             (1+c[2]*exp(c[3]*(x[i]-xc)))
             -e[0]*exp(e[1]*x[i])*c[0]*exp(c[1]*(-x[i]+xb))*(-x[i]+xb)/
             pow(1+c[0]*exp(c[1]*(-x[i]+xb)), 2);
    df_dc3 = -gsl_fit_poly_eval(x[i], p, order)*
             c[2]*exp(c[3]*(x[i]-xc))*(x[i]-xc)/
             (1+c[0]*exp(c[1]*(x[i]-xb)))/
             pow(1+c[2]*exp(c[3]*(x[i]-xc)), 2)
             -(E[0]*exp(E[1]*x[i])+E[2])*c[2]*exp(c[3]*(-x[i]+xc))*(-x[i]+xc)/
             pow(1+c[2]*exp(c[3]*(-x[i]+xc)), 2);
    gsl_matrix_set(J, i, 0, df_dc0);
    gsl_matrix_set(J, i, 1, df_dc1);
    gsl_matrix_set(J, i, 2, df_dc2);
    gsl_matrix_set(J, i, 3, df_dc3);
  }

  return GSL_SUCCESS;
}

int gsl_fit_epe_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                    gsl_matrix* J)
{
  gsl_fit_epe_f(coeff, params, f);
  gsl_fit_epe_df(coeff, params, J);

  return GSL_SUCCESS;
}

double gsl_fit_epe_eval_f(double x, double xc, double xb, double* c, double* p,
                          int order, double* e, double* E)
{

  return gsl_fit_poly_eval(x, p, order)/
         (1+c[0]*exp(c[1]*(x-xb)))/
         (1+c[2]*exp(c[3]*(x-xc)))+
         e[0]*exp(e[1]*x)/(1+c[0]*exp(c[1]*(-x+xb)))+
         (E[0]*exp(E[1]*x)+E[2])/(1+c[2]*exp(c[3]*(-x+xc)));
}

double gsl_fit_epe_eval_df(double x, double xc, double xb, double* c,
                           double* p, double *dp, int order,
                           double* e, double* E)
{
  double nom1, nom2, nom3, den1, den2, den3;

  nom1 = (gsl_fit_poly_eval(x, dp, order-1)*(1+c[0]*exp(c[1]*(x-xb)))*
         (1+c[2]*exp(c[3]*(x-xc)))-
         gsl_fit_poly_eval(x, p, order)*((1+c[0]*exp(c[1]*(x-xb)))*
         c[2]*c[3]*exp(c[3]*(x-xc))+
         c[0]*c[1]*exp(c[1]*(x-xb))*(1+c[2]*exp(c[3]*(x-xc)))));
  den1 = pow(1+c[0]*exp(c[1]*(x-xb)), 2)*pow(1+c[2]*exp(c[3]*(x-xc)), 2);

  nom2 = (e[0]*e[1]*exp(e[1]*x)*(1+c[0]*exp(c[1]*(-x+xb)))+
         e[0]*exp(e[1]*x)*c[0]*c[1]*exp(c[1]*(-x+xb)));
  den2 = pow(1+c[0]*exp(c[1]*(-x+xb)), 2);

  nom3 = (E[0]*E[1]*exp(E[1]*x)*(1+c[2]*exp(c[3]*(-x+xc)))+
         (E[0]*exp(E[1]*x)+E[2])*c[2]*c[3]*exp(c[3]*(-x+xc)));
  den3 = pow(1+c[2]*exp(c[3]*(-x+xc)), 2);

  double res = 0;

  if (den1 != numeric_limits<double>::infinity()) res += nom1/den1;
  if (den2 != numeric_limits<double>::infinity()) res += nom2/den2;
  if (den3 != numeric_limits<double>::infinity()) res += nom3/den3;

  return res;
}


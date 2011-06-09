
#ifndef GSR_FIT
#define GSR_FIT

// library headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

struct gsr_fit_data {
  double *x, *y;
  double *p, *e, *E;
  int order;
  double xc, xb;
  double *sigma;
  int n;
};

int gsr_fit_f(const gsl_vector* coeff, void* params, gsl_vector* f)
{
  double *x = ((gsr_fit_data*)params)->x;
  double *y = ((gsr_fit_data*)params)->y;
  double *p = ((gsr_fit_data*)params)->p;
  double *e = ((gsr_fit_data*)params)->e;
  double *E = ((gsr_fit_data*)params)->E;
  int order = ((gsr_fit_data*)params)->order;
  double xc = ((gsr_fit_data*)params)->xc;
  double xb = ((gsr_fit_data*)params)->xb;
  double *sigma = ((gsr_fit_data*)params)->sigma;
  int n = ((gsr_fit_data*)params)->n;

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
    gsl_vector_set(f, i, (Yi-y[i])/sigma[i]);
  }

  return GSL_SUCCESS;
}

int gsr_fit_df(const gsl_vector* coeff, void* params, gsl_matrix* J)
{
  double *x = ((gsr_fit_data*)params)->x;
  double *y = ((gsr_fit_data*)params)->y;
  double *p = ((gsr_fit_data*)params)->p;
  double *e = ((gsr_fit_data*)params)->e;
  double *E = ((gsr_fit_data*)params)->E;
  int order = ((gsr_fit_data*)params)->order;
  double xc = ((gsr_fit_data*)params)->xc;
  double xb = ((gsr_fit_data*)params)->xb;
  double *sigma = ((gsr_fit_data*)params)->sigma;
  int n = ((gsr_fit_data*)params)->n;

  const int nc = 4;
  double c[nc];
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double df_dc1, df_dc2, df_dc3, df_dc4;
  for (int i = 0; i < n; i++) {
    df_dc1 = -gsl_fit_poly_eval(x[i], p, order)*exp(c[1]*(x[i]-xb))/
             pow(1+c[0]*exp(c[1]*(x[i]-xb)), 2)/
             (1+c[2]*exp(c[3]*(x[i]-xc)))
             -e[0]*exp(e[1]*x[i])*exp(c[1]*(-x[i]+xb))/
             pow(1+c[0]*exp(c[1]*(-x[i]+xb)), 2);
    df_dc2 = -gsl_fit_poly_eval(x[i], p, order)*exp(c[3]*(x[i]-xc))/
             (1+c[0]*exp(c[1]*(x[i]-xb)))/
             pow(1+c[2]*exp(c[3]*(x[i]-xc)), 2)
             -(E[0]*exp(E[1]*x[i])+E[2])*exp(c[3]*(-x[i]+xc))/
             pow(1+c[2]*exp(c[3]*(-x[i]+xc)), 2);
    df_dc3 = -gsl_fit_poly_eval(x[i], p, order)*
             c[0]*exp(c[1]*(x[i]-xb))*(x[i]-xb)/
             pow(1+c[0]*exp(c[1]*(x[i]-xb)), 2)/
             (1+c[2]*exp(c[3]*(x[i]-xc)))
             -e[0]*exp(e[1]*x[i])*c[0]*exp(c[1]*(-x[i]+xb))*(-x[i]+xb)/
             pow(1+c[0]*exp(c[1]*(-x[i]+xb)), 2);
    df_dc4 = -gsl_fit_poly_eval(x[i], p, order)*
             c[2]*exp(c[3]*(x[i]-xc))*(x[i]-xc)/
             (1+c[0]*exp(c[1]*(x[i]-xb)))/
             pow(1+c[2]*exp(c[3]*(x[i]-xc)), 2)
             -(E[0]*exp(E[1]*x[i])+E[2])*c[2]*exp(c[3]*(-x[i]+xc))*(-x[i]+xc)/
             pow(1+c[2]*exp(c[3]*(-x[i]+xc)), 2);
    gsl_matrix_set(J, i, 0, df_dc1/sigma[i]);
    gsl_matrix_set(J, i, 1, df_dc2/sigma[i]);
    gsl_matrix_set(J, i, 2, df_dc3/sigma[i]);
    gsl_matrix_set(J, i, 3, df_dc4/sigma[i]);
  }

  return GSL_SUCCESS;
}

int gsr_fit_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                gsl_matrix* J)
{
  gsr_fit_f(coeff, params, f);
  gsr_fit_df(coeff, params, J);

  return GSL_SUCCESS;
}

double gsr_fit_eval_f(double x, double xc, double xb, double* c, double* p,
                      int order, double* e, double* E)
{

  return gsl_fit_poly_eval(x, p, order)/
         (1+c[0]*exp(c[1]*(x-xb)))/
         (1+c[2]*exp(c[3]*(x-xc)))+
         e[0]*exp(e[1]*x)/(1+c[0]*exp(c[1]*(-x+xb)))+
         (E[0]*exp(E[1]*x)+E[2])/(1+c[2]*exp(c[3]*(-x+xc)));
}

double gsr_fit_eval_df(double x, double xc, double xb, double* c,
                       double* p, double *dp, int order,
                       double* e, double* E)
{

  return (gsl_fit_poly_eval(x, dp, order-1)*(1+c[0]*exp(c[1]*(x-xb)))*
         (1+c[2]*exp(c[3]*(x-xc)))-
         gsl_fit_poly_eval(x, p, order)*((1+c[0]*exp(c[1]*(x-xb)))*
         c[2]*c[3]*exp(c[3]*(x-xc))+
         c[0]*c[1]*exp(c[1]*(x-xb))*(1+c[2]*exp(c[3]*(x-xc)))))/
         pow(1+c[0]*exp(c[1]*(x-xb)), 2)/pow(1+c[2]*exp(c[3]*(x-xc)), 2)+
         (e[0]*e[1]*exp(e[1]*x)*(1+c[0]*exp(c[1]*(-x+xb)))+
         e[0]*exp(e[1]*x)*c[0]*c[1]*exp(c[1]*(-x+xb)))/
         pow(1+c[0]*exp(c[1]*(-x+xb)), 2)+
         (E[0]*E[1]*exp(E[1]*x)*(1+c[2]*exp(c[3]*(-x+xc)))+
         (E[0]*exp(E[1]*x)+E[2])*c[2]*c[3]*exp(c[3]*(-x+xc)))/
         pow(1+c[2]*exp(c[3]*(-x+xc)), 2);
}

#endif


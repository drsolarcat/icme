
// project headers
#include "gsl_fit_exp.h"
#include "curve.h"
#include "integrator.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
// standard headers
#include <iostream>

using namespace std;
using namespace Eigen;

int gsl_fit_exp_0(int n, double* x, double* y, double* c0)
{
  Map<VectorXd> xVec(x, n);
  Map<VectorXd> yVec(y, n);

  double xMin = xVec.minCoeff(), xMax = xVec.maxCoeff();
  double X = xMax-xMin;
  int m = 200;

  Curve curve(xVec, yVec);
  curve.resample(m);

  Integrator integrator;

  double I1, I2, I3, I4;

  I1 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(0)*curve.cols().y.array())/1;
  I2 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(1)*curve.cols().y.array())/1;
  I3 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(2)*curve.cols().y.array())/2;
  I4 = integrator.NewtonCotes(2, curve.cols().x,
    (xMax-curve.cols().x.array()).pow(3)*curve.cols().y.array())/6;

  double tau = (12*I4-6*I3*X+I2*pow(X, 2))/(-12*I3+6*I2*X-I1*pow(X, 2));
  double Q1 = exp(-xMin/tau);
  double Q = exp(-X/tau);
  c0[0] = 2/X/((1+Q)*X+2*(Q-1)*tau)*(I2*(Q-1)+I1*(X+(Q-1)*tau));
  c0[1] = (2*I2-I1*X)/tau/((1+Q)*X+2*(Q-1)*tau)/Q1;
  c0[2] = -1/tau;

  return GSL_SUCCESS;
}

int gsl_fit_exp(int n, double* x, double* y, double* c)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  const int nc = 3;

  double c0[3];
  gsl_fit_exp_0(n, x, y, c0);

  gsl_multifit_function_fdf f;
  gsl_fit_exp_data params;

  params.x = x;
  params.y = y;
  params.n = n;

  f.f = &gsl_fit_exp_f;
  f.df = &gsl_fit_exp_df;
  f.fdf = &gsl_fit_exp_fdf;
  f.n = n;
  f.p = nc;
  f.params = &params;

  gsl_vector *coeff0 = gsl_vector_alloc(nc);
  gsl_vector_set(coeff0, 0, c0[0]);
  gsl_vector_set(coeff0, 1, c0[1]);
  gsl_vector_set(coeff0, 2, c0[2]);

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, nc);
  gsl_multifit_fdfsolver_set(s, &f, coeff0);

  int i = 0;
  do {
    i++;
    status = gsl_multifit_fdfsolver_iterate(s);
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-40, 1e-40);
  } while (status == GSL_CONTINUE && i < 5000);

  c[0] = gsl_vector_get(s->x, 0);
  c[1] = gsl_vector_get(s->x, 1);
  c[2] = gsl_vector_get(s->x, 2);

  gsl_multifit_fdfsolver_free (s);

  return GSL_SUCCESS;
}

int gsl_fit_exp_f(const gsl_vector* coeff, void* params, gsl_vector* f)
{
  double *x = ((gsl_fit_exp_data*)params)->x;
  double *y = ((gsl_fit_exp_data*)params)->y;
  int n = ((gsl_fit_exp_data*)params)->n;

  const int nc = 3;
  double c[nc];
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double Yi;

  for (int i = 0; i < n; i++) {
    Yi = c[0]+c[1]*exp(c[2]*x[i]);
    gsl_vector_set(f, i, Yi-y[i]);
  }

  return GSL_SUCCESS;
}

int gsl_fit_exp_df(const gsl_vector* coeff, void* params, gsl_matrix* J)
{
  double *x = ((gsl_fit_exp_data*)params)->x;
  double *y = ((gsl_fit_exp_data*)params)->y;
  int n = ((gsl_fit_exp_data*)params)->n;

  const int nc = 3;
  double c[nc];
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double df_dc0, df_dc1, df_dc2;
  for (int i = 0; i < n; i++) {
    df_dc0 = 1;
    df_dc1 = exp(c[2]*x[i]);
    df_dc2 = c[1]*x[i]*exp(c[2]*x[i]);

    gsl_matrix_set(J, i, 0, df_dc0);
    gsl_matrix_set(J, i, 1, df_dc1);
    gsl_matrix_set(J, i, 2, df_dc2);
  }

  return GSL_SUCCESS;
}

int gsl_fit_exp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                    gsl_matrix* J)
{
  gsl_fit_exp_f(coeff, params, f);
  gsl_fit_exp_df(coeff, params, J);

  return GSL_SUCCESS;
}

double gsl_fit_exp_eval_f(double x, double* c)
{
  return c[0]+c[1]*exp(c[2]*x);
}

double gsl_fit_exp_eval_df(double x, double* c)
{
  return c[1]*c[2]*exp(c[2]*x);
}


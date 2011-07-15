
// project headers
#include "fit_exp.h"
#include "fit_cexp.h"
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

// process the fitting
int fit_exp(int n, double* x, double* y, double* c)
{
  const gsl_multifit_fdfsolver_type *T; // solver type
  gsl_multifit_fdfsolver *s; // solver
  int status;
  const int nc = 2; // number of fitting coefficients

  double c0[2]; // initialize initial guess for fitting parameters
  fit_exp_0(n, x, y, c0); // compute the initial guess

  gsl_multifit_function_fdf f; // fitting function
  fit_exp_data params; // fitting function parameters

  // initialize parametersof the fitting function
  params.x = x;
  params.y = y;
  params.n = n;

  // initialize the fitting function
  f.f = &fit_exp_f; // function
  f.df = &fit_exp_df; // Jacobian
  f.fdf = &fit_exp_fdf; // function and Jacobian
  f.n = n; // length
  f.p = nc; // number of fitting coefficients
  f.params = &params; // fitting function parameters

  // vector for initial guess coefficients
  gsl_vector* cVec0 = gsl_vector_alloc(nc);
  gsl_vector_set(cVec0, 0, c0[0]);
  gsl_vector_set(cVec0, 1, c0[1]);

  // set the solver
  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, nc);
  gsl_multifit_fdfsolver_set(s, &f, cVec0);

  // iterate the solver
  int i = 0;
  do {
    i++;
    status = gsl_multifit_fdfsolver_iterate(s);
    status = gsl_multifit_test_delta(s->dx, s->x, 1e-40, 1e-40);
  } while (status == GSL_CONTINUE && i < 5000);

  // coefficients of the fitted curve
  c[0] = gsl_vector_get(s->x, 0);
  c[1] = gsl_vector_get(s->x, 1);

  // free the solver
  gsl_multifit_fdfsolver_free (s);

  return GSL_SUCCESS; // done
}

// compute the initial guess for the fitting parameters
int fit_exp_0(int n, double* x, double* y, double* c0)
{
  double cc0[3];
  fit_cexp_0(n, x, y, cc0);
  c0[0] = cc0[1];
  c0[1] = cc0[2];

  return GSL_SUCCESS; // done
}

// minimization function
int fit_exp_f(const gsl_vector* coeff, void* params, gsl_vector* f)
{
  // get function parameters
  double* x = ((fit_exp_data*)params)->x;
  double* y = ((fit_exp_data*)params)->y;
  int n = ((fit_exp_data*)params)->n;

  const int nc = 2; // number of fitting coefficients
  double c[nc]; // coefficients array
  // fill the coefficients array
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double Yi; // fitting function value

  // fill the minimization function values
  for (int i = 0; i < n; i++) {
    Yi = c[0]*exp(c[1]*x[i]);
    gsl_vector_set(f, i, Yi-y[i]);
  }

  return GSL_SUCCESS; // done
}

// compute Jacobian of the fitting function
int fit_exp_df(const gsl_vector* coeff, void* params, gsl_matrix* J)
{
  // get function parameters from the parameters structure
  double* x = ((fit_exp_data*)params)->x;
  double* y = ((fit_exp_data*)params)->y;
  int n = ((fit_exp_data*)params)->n;

  const int nc = 2; // number of  fitting coefficients
  double c[nc]; // initialize array of the fitting coefficients
  // fill the coefficients array
  for (int i = 0; i < nc; i++) {
    c[i] = gsl_vector_get(coeff, i);
  }

  double df_dc0, df_dc1, df_dc2; // derivatives over fitting coefficients
  // fill the Jacobian with derivatives
  for (int i = 0; i < n; i++) {
    df_dc0 = exp(c[1]*x[i]); // df/dc0
    df_dc1 = c[0]*x[i]*exp(c[1]*x[i]); // df/dc1

    gsl_matrix_set(J, i, 0, df_dc0);
    gsl_matrix_set(J, i, 1, df_dc1);
  }

  return GSL_SUCCESS; // done
}

// compute Jacobian and function simultaiously
int fit_exp_fdf(const gsl_vector* coeff, void* params, gsl_vector* f,
                gsl_matrix* J)
{
  fit_exp_f(coeff, params, f); // compute the minimization function
  fit_exp_df(coeff, params, J); // compute the Jacobian

  return GSL_SUCCESS; // done
}

// evaluate the fitting function
double fit_exp_eval_f(double x, double* c)
{
  return c[0]*exp(c[1]*x); // return the f(x) value
}

// evaluate the derivative of the fitting function
double fit_exp_eval_df(double x, double* c)
{
  return c[0]*c[1]*exp(c[1]*x); // return the df/dx(x) value
}


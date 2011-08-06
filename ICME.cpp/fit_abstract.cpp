
// project headers
#include "fit_abstract.h"
// library headers
#include <gsl/gsl_multifit_nlin.h>
// standard headers
#include <iostream>

using namespace std;

// define functions for fitting
extern "C"
{
  // calculates residues
  int fit_f(const gsl_vector* coeff, void* params, gsl_vector* f) {
    AbstractFit *fit = (AbstractFit*)params;
    return fit->f(coeff, params, f);
  }
  // calculates Jacobian
  int fit_df(const gsl_vector* coeff, void* params, gsl_matrix* J) {
    AbstractFit *fit = (AbstractFit*)params;
    return fit->df(coeff, params, J);
  }
  // calculates residues and Jacobian
  int fit_fdf(const gsl_vector* coeff, void* params,
              gsl_vector* f, gsl_matrix* J)
  {
    AbstractFit *fit = (AbstractFit*)params;
    return fit->fdf(coeff, params, f, J);
  }
}

// fit method
void AbstractFit::fit()
{
  // status flag
  int status;

  // compute the initial guess
  fit0();

  gsl_multifit_function_fdf f; // fitting function

  // initialize the fitting function
  f.f = &fit_f; // function
  f.df = &fit_df; // Jacobian
  f.fdf = &fit_fdf; // function and Jacobian
  f.n = _n; // length
  f.p = _nc; // number of fitting coefficients
  f.params = this; // fitting function parameters

  // vector for initial guess coefficients
  gsl_vector *cVec0 = gsl_vector_alloc(_nc);
  for (int i = 0; i < _nc; i++) {
    gsl_vector_set(cVec0, i, _c0[i]);
  }

  // initialize the solver
  gsl_multifit_fdfsolver *solver =
    gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, _n, _nc);
  gsl_multifit_fdfsolver_set(solver, &f, cVec0);

  // iterate the solver
  int i = 0;
  do {
    i++;
    status = gsl_multifit_fdfsolver_iterate(solver);
    status = gsl_multifit_test_delta(solver->dx, solver->x, _epsabs, _epsrel);
  } while (status == GSL_CONTINUE && i < _niter);

  // coefficients of the fitted curve
  for (int i = 0; i < _nc; i++) {
    _c[i] = gsl_vector_get(solver->x, i);
  }

  // free the solver
  gsl_multifit_fdfsolver_free (solver);
}

// calculate vector of residues between fitted function and data
int AbstractFit::f(const gsl_vector* cVec, void* params, gsl_vector* fVec)
{
  // get function parameters
  AbstractFit *fit = (AbstractFit*)params;

  double c[fit->nc()]; // coefficients array
  // fill the coefficients array
  for (int i = 0; i < fit->nc(); i++) {
    c[i] = gsl_vector_get(cVec, i);
  }

  double Yi; // fitting function value
  // fill the minimization function values
  for (int i = 0; i < fit->n(); i++) {
    gsl_vector_set(fVec, i, fit->f(fit->x()[i], c)-fit->y()[i]);
  }

  return GSL_SUCCESS; // done
}

// calculate Jacobian
int AbstractFit::df(const gsl_vector* cVec, void* params, gsl_matrix* J)
{
  // get function parameters from the parameters structure
  AbstractFit *fit = (AbstractFit*)params;

  double c[fit->nc()]; // initialize array of the fitting coefficients
  // fill the coefficients array
  for (int i = 0; i < fit->nc(); i++) {
    c[i] = gsl_vector_get(cVec, i);
  }

  double dy[fit->nc()]; // derivatives over fitting coefficients
  // fill the Jacobian with derivatives
  for (int i = 0; i < fit->n(); i++) {
    df(fit->x()[i], c, dy);
    for (int k = 0; k < fit->nc(); k++) {
      gsl_matrix_set(J, i, k, dy[k]);
    }
  }

  return GSL_SUCCESS; // done
}

// calculate both residues and Jacobian
int AbstractFit::fdf(const gsl_vector* cVec, void* params,
                     gsl_vector* fVec, gsl_matrix* J)
{
  f(cVec, params, fVec); // compute the minimization function
  df(cVec, params, J); // compute the Jacobian

  return GSL_SUCCESS; // done
}

// evaluate function at all points in array X
void AbstractFit::f(int n, double* x, double* y)
{
  for (int i = 0; i < n; i++) {
    y[i] = f(x[i]);
  }
}

// evaluate function derivative at all points in array X
void AbstractFit::df(int n, double* x, double* dy)
{
  for (int i = 0; i < n; i++) {
    dy[i] = df(x[i]);
  }
}


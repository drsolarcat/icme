
#ifndef MYBLAS_H
#define MYBLAS_H

using namespace std;

int gsl_blas_dcross(const gsl_vector *, const gsl_vector *, gsl_vector *);
int gsl_blas_drot(const gsl_vector *, const gsl_vector *, double, gsl_vector *);
gsl_vector * vcross(const gsl_vector *, const gsl_vector *);
double vdot(const gsl_vector * x, const gsl_vector * y);
gsl_vector * vrot(const gsl_vector * x, const gsl_vector * a, double phi);
gsl_vector * vadd(const gsl_vector * x, const gsl_vector * y);
gsl_vector * vsub(const gsl_vector * x, const gsl_vector * y);
gsl_vector * vscale(const gsl_vector *, const double);
gsl_vector * vnormalize(const gsl_vector *);


int gsl_blas_dcross(
  const gsl_vector * x,
  const gsl_vector * y,
  gsl_vector * result)
{
  gsl_vector_set(result, 0,
    gsl_vector_get(x, 2)*gsl_vector_get(y, 3)-
    gsl_vector_get(x, 3)*gsl_vector_get(y, 2));
  gsl_vector_set(result, 1,
    gsl_vector_get(x, 3)*gsl_vector_get(y, 1)-
    gsl_vector_get(x, 1)*gsl_vector_get(y, 3));
  gsl_vector_set(result, 2,
    gsl_vector_get(x, 1)*gsl_vector_get(y, 2)-
    gsl_vector_get(x, 2)*gsl_vector_get(y, 1));
  return GSL_SUCCESS;
}

int gsl_blas_drot(
  const gsl_vector * x,
  const gsl_vector * a,
  double phi,
  gsl_vector * result)
{
  gsl_vector * x_dot_a;
  gsl_blas_ddot(x, a, xdota);
  gsl_blas_dcross(x, a, xcrossa);
  for(int i=0; i<3; i++)
    gsl_vector_set(result, i,
      gsl_vector_get(x, i)*cos(phi)+
      xdota*(1-cos(phi))*gsl_vector_get(a, i)-
      gsl_vector_get(xcrossa, i)*sin(phi));
  return GSL_SUCCESS;
}

gsl_vector * vcross(const gsl_vector * x, const gsl_vector * y) {
  gsl_vector * result;
  gsl_blas_dcross(x, y, result);
  return result;
}

double vdot(const gsl_vector * x, const gsl_vector * y) {
  double result;
  gsl_blas_ddot(x, y, &result);
  return result;
}

gsl_vector * vrot(const gsl_vector * x, const gsl_vector * a, double phi) {
  gsl_vector * result;
  gsl_blas_drot(x, a, phi, result);
  return result;
}

gsl_vector * vadd(const gsl_vector * x, const gsl_vector * y) {
  gsl_vector * result;
  gsl_vector_memcpy(result, x);
  gsl_vector_add(result, y);
  return result;
}

gsl_vector * vsub(const gsl_vector * x, const gsl_vector * y) {
  gsl_vector * result;
  gsl_vector_memcpy(result, x);
  gsl_vector_sub(result, y);
  return result;
}

gsl_vector * vscale(const gsl_vector * x, const double a) {
  gsl_vector * result;
  gsl_vector_memcpy(result, x);
  gsl_vector_scale(result, a);
  return result;
}

gsl_vector * vnormalize(const gsl_vector * x) {
  double norm = gsl_blas_dnrm2(x);
  return vscale(x, 1/norm);
}

#endif


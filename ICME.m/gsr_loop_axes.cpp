
#include <vector>
#include <gsl/gsl_const_mksa.h>
#include <cmath>
#include "myblas.h"
#include "mex.h"

using namespace std;

double trapezoid_integrate(vector<double> x, vector<double> y) {
  double sum = 0;

  for(int i=1; i<x.size(); i++)
    sum += (x[i]-x[i-1])*(y[i]-y[i-1]);
  return sum;
}

void projection(vector<double> &Ax, vector<double> &Ay, vector<double> &Az,
  vector<double> x, vector<double> y, vector<double> z)
{
  vector<double> A(3);

  for(int i=0; i<Ax.size(); i++) {
    A[0] = Ax[i];
    A[1] = Ay[i];
    A[2] = Az[i];
    Ax[i] = vdot(A, x);
    Ay[i] = vdot(A, y);
    Az[i] = vdot(A, z);
  }
}

void computeCurve(vector<double> Vht,
  vector<double> Bx, vector<double> By, vector<double> Bz, vector<double> Pth,
  double dt, vector<double> x, vector<double> y, vector<double> z)
{
  int n = Bx.size();
  vector<double> A(n), Pt(n), X(n), Y(n);
  double dx = -vdot(Vht, x)*dt;

  projection(Bx, By, Bz, x, y, z);

  for(int i=0; i<n; i++) {
    X.clear();
    Y.clear();
    for(int k=0; k<i; k++) {
      X.push_back(dx*k);
      Y.push_back(By[k]);
    }
    A.push_back(i == 0 ? 0 : trapezoid_integrate(X, Y));
    Pt.push_back(Pth[i]+pow(Bz[i], 2)/2/GSL_CONST_MKSA_VACUUM_PERMEABILITY);
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double minTheta, dTheta, maxTheta, minPhi, dPhi, maxPhi, theta, phi, dt;
  int n;
  vector<double> z1(3), x2(3), y2(3), z2(3);

  minTheta = mxGetScalar(prhs[0]);
  dTheta   = mxGetScalar(prhs[1]);
  maxTheta = mxGetScalar(prhs[2]);
  minPhi   = mxGetScalar(prhs[3]);
  dPhi     = mxGetScalar(prhs[4]);
  maxPhi   = mxGetScalar(prhs[5]);
  vector<double> x(mxGetPr(prhs[6]), mxGetPr(prhs[6])+sizeof(double)*3);
  vector<double> y(mxGetPr(prhs[7]), mxGetPr(prhs[7])+sizeof(double)*3);
  vector<double> z(mxGetPr(prhs[8]), mxGetPr(prhs[8])+sizeof(double)*3);
  double aVht [] = {mxGetScalar(prhs[9]), mxGetScalar(prhs[10]),
    mxGetScalar(prhs[11])};
  vector<double> Vht(aVht, aVht+sizeof(double)*3);
  n = mxGetN(prhs[12]);
  vector<double>  Bx(mxGetPr(prhs[12]), mxGetPr(prhs[12])+sizeof(double)*n);
  vector<double>  By(mxGetPr(prhs[13]), mxGetPr(prhs[13])+sizeof(double)*n);
  vector<double>  Bz(mxGetPr(prhs[14]), mxGetPr(prhs[14])+sizeof(double)*n);
  vector<double> Pth(mxGetPr(prhs[15]), mxGetPr(prhs[15])+sizeof(double)*n);
  dt = mxGetScalar(prhs[16]);

  int nTheta = round((maxTheta-minTheta)/dTheta+1);
  int nPhi = round((maxPhi-minPhi)/dPhi+1);

  double residual[nTheta][nPhi];

  theta = minTheta;
  while(theta <= maxTheta) {
    z1 = vrotate(z, y, theta*M_PI/180);
    phi = minPhi;
    while(phi <= maxPhi) {
      z2 = vrotate(z1, z, phi*M_PI/180);
      x2 = vnormalize(-(Vht-vdot(Vht, z2)*z2));
      y2 = vcross(z2, x2);
      computeCurve(Vht, Bx, By, Bz, Pth, dt, x2, y2, z2);
      /*
      now we should get the curve and calculate its residue
      possible solution:
      1. get the branches in C++
      2. calculate residue in C++
      3. everything else is done in Matlab
      */

      phi += dPhi;
    }
    theta += dTheta;
  }
}


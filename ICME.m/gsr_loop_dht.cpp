
#include <cmath>
#include <cstring>
#include <vector>
#include "myblas.h"
#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  double minVx, dVx, maxVx, minVy, dVy, maxVy, minVz, dVz, maxVz;
  double Vfx, Vfy, Vfz;
  double Vhtx, Vhty, Vhtz;
  double D = 0, minD = 0;
  double *Vx, *Vy, *Vz, *Bx, *By, *Bz;
  int n;
  vector<double> a(3), b(3), c(3);

  minVx = mxGetScalar(prhs[0]);
  dVx   = mxGetScalar(prhs[1]);
  maxVx = mxGetScalar(prhs[2]);
  minVy = mxGetScalar(prhs[3]);
  dVy   = mxGetScalar(prhs[4]);
  maxVy = mxGetScalar(prhs[5]);
  minVz = mxGetScalar(prhs[6]);
  dVz   = mxGetScalar(prhs[7]);
  maxVz = mxGetScalar(prhs[8]);
  Vx    = mxGetPr(prhs[9]);
  Vy    = mxGetPr(prhs[10]);
  Vz    = mxGetPr(prhs[11]);
  Bx    = mxGetPr(prhs[12]);
  By    = mxGetPr(prhs[13]);
  Bz    = mxGetPr(prhs[14]);
  n     = mxGetN(prhs[9]);

  Vhtx = minVx;
  Vhty = minVy;
  Vhtz = minVz;

  Vfx = minVx;
  while(Vfx <= maxVx) {
    Vfy = minVy;
    while(Vfy <= maxVy) {
      Vfz = minVz;
      while(Vfz <= maxVz) {
        for(int i=0; i<n; i++) {
          a[0] = Vx[i]-Vfx;
          a[1] = Vy[i]-Vfy;
          a[2] = Vz[i]-Vfz;
          b[0] = Bx[i];
          b[1] = By[i];
          b[2] = Bz[i];
          c = vcross(a, b);
          D += pow(c[0], 2)+pow(c[1], 2)+pow(c[2], 2);
        }
        D /= n;
        if (minD == 0) {
          minD = D;
          continue;
        }
        if (D < minD) {
          minD = D;
          Vhtx = Vfx;
          Vhty = Vfy;
          Vhtz = Vfz;
        }
        Vfz += dVz;
      }
      Vfy += dVy;
    }
    Vfx += dVx;
  }

  double *pointers[] = {&Vhtx, &Vhty, &Vhtz};
  for(int i=0; i<=2; i++) {
    plhs[i] = mxCreateDoubleMatrix(1, 1, mxREAL);
    memcpy(mxGetPr(plhs[i]), pointers[i], sizeof(double));
  }
}

/*
decide wether to write everything in C++

functions we will need (simple(1) - complicated(5)):
2 - cross
2 - dot
4 - corrcoef
2 - norm
5 - eig
5 - smooth
5 - lsqcurvefit
4 - polyfit
3 - polyder
4 - sortrows
5 - resample
*/


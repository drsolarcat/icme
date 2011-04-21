
#include "mex.h"
#include <vector>
#include <cstring>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int n = 1000000;
  vector<double> outVector;

  for (int i=0; i<n; i++)
    outVector.push_back(double(i));

  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  for (int i=0; i<10; i++)
    memcpy(mxGetPr(plhs[0]), &outVector[0], sizeof(double)*outVector.size());
}


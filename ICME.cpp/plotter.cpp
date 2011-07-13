
// standard headers
#include <iostream>
// project headers
#include "plotter.h"

using namespace std;
using namespace Eigen;

Plotter::Plotter() {
  _matlab = engOpen(NULL);
  engOutputBuffer(_matlab, NULL, 0);
}

Plotter::~Plotter() {
//  engClose(_matlab);
}

void Plotter::plotResidueMap(const MatrixXd& residue, const VectorXd& theta,
                             const VectorXd& phi,
                             double optTheta, double optPhi) {

  mxArray *mResidue, *mTheta, *mPhi, *mOptTheta, *mOptPhi;

  mTheta = mxCreateDoubleMatrix(1, theta.size(), mxREAL);
  mPhi = mxCreateDoubleMatrix(1, phi.size(), mxREAL);

  mOptTheta = mxCreateDoubleScalar(optTheta);
  mOptPhi = mxCreateDoubleScalar(optPhi);

  memcpy((double*)mxGetPr(mTheta), (double*)theta.data(),
         theta.size()*sizeof(double));
  memcpy((double*)mxGetPr(mPhi), (double*)phi.data(),
         phi.size()*sizeof(double));

  engPutVariable(_matlab, "theta", mTheta);
  engPutVariable(_matlab, "phi", mPhi);
  engPutVariable(_matlab, "optTheta", mOptTheta);
  engPutVariable(_matlab, "optPhi", mOptPhi);

  mResidue = mxCreateDoubleMatrix(residue.rows(), residue.cols(), mxREAL);

  memcpy((double*)mxGetPr(mResidue), (double*)residue.data(),
         residue.rows()*residue.cols()*sizeof(double));

  engPutVariable(_matlab, "R", mResidue);

  char matlabCommand[] =
  "R = R.^-1;\n"
  "figure\n"
  "h = polar([0 2*pi], [0 90]);\n"
  "ph = findall(gca, 'type', 'patch');\n"
  "set(ph, 'facecolor', [0 0 .5], 'edgecolor', [0 0 .5]);\n"
  "set(gcf, 'Color', 'white');\n"
  "set(gcf, 'Renderer', 'Painters');\n"
  "pl = findobj(allchild(gca));\n"
  "hold on\n"
  "contour(theta'*cos(phi*pi/180), theta'*sin(phi*pi/180), R, ...\n"
  "        linspace(min(R(:)), max(R(:)), 500), 'Fill', 'on');\n"
  "plot(optTheta*cos(optPhi*pi/180), optTheta*sin(optPhi*pi/180), ...\n"
  "  '.w', 'MarkerSize', 20);\n"
  "colorbar\n"
  "for i = 1:length(pl)-1\n"
  "  if strcmpi(get(pl(i), 'Type'), 'line')\n"
  "    set(pl(i), 'Color', 'white');\n"
  "  elseif strcmpi(get(pl(i), 'Type'), 'text') && i > 25\n"
  "    set(pl(i), 'Color', 'white');\n"
  "  end\n"
  "  uistack(pl(i), 'top');\n"
  "end\n"
  "delete(h)\n"
  "cbar_handle = findobj(gcf, 'Tag', 'Colorbar');\n"
  "set(get(cbar_handle,'xlabel'),'string','1/R');\n"
  "hold off\n";

  engEvalString(_matlab, matlabCommand);
}

void Plotter::plotMagneticMap(const MatrixXd& Axy, const MatrixXd& Bz,
                              const VectorXd& X, const VectorXd& Y) {

  mxArray *mAxy, *mBz, *mX, *mY;
  mX = mxCreateDoubleMatrix(1, X.size(), mxREAL);
  mY = mxCreateDoubleMatrix(1, Y.size(), mxREAL);

  memcpy((double*)mxGetPr(mX), (double*)X.data(), X.size()*sizeof(double));
  memcpy((double*)mxGetPr(mY), (double*)Y.data(), Y.size()*sizeof(double));

cout << Bz.minCoeff() << endl << Bz.maxCoeff() << endl;
  engPutVariable(_matlab, "X", mX);
  engPutVariable(_matlab, "Y", mY);

  mAxy = mxCreateDoubleMatrix(Axy.rows(), Axy.cols(), mxREAL);

  memcpy((double*)mxGetPr(mAxy), (double*)Axy.data(),
         Axy.rows()*Axy.cols()*sizeof(double));

  engPutVariable(_matlab, "Axy", mAxy);

  mBz = mxCreateDoubleMatrix(Bz.rows(), Bz.cols(), mxREAL);

  memcpy((double*)mxGetPr(mBz), (double*)Bz.data(),
         Bz.rows()*Bz.cols()*sizeof(double));

  engPutVariable(_matlab, "Bz", mBz);

  char matlabCommand[] =
  "figure\n"
  "contour(X, Y, Bz, linspace(min(Bz(:)), max(Bz(:)), 100), ...\n"
  "  'Fill', 'on', 'LineStyle', 'none');\n"
  "colorbar\n"
  "hold on\n"
  "contour(X, Y, Axy, 50, 'LineColor', 'black');\n"
  "caxis([min(Bz(:)) max(Bz(:))]);\n"
  "set(gca, 'Color', [0 0 .5]);\n"
  "set(gcf, 'Color', 'white');\n"
  "set(gcf, 'InvertHardCopy', 'off');\n"
  "xlabel 'X_{MC} (AU)'\n"
  "ylabel 'Y_{MC} (AU)'\n";

  engEvalString(_matlab, matlabCommand);
}


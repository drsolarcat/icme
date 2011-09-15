
// project headers
#include "plotter.h"
#include "gnuplot.h"
// standard headers
#include <iostream>
#include <string>
#include <sys/stat.h>

using namespace std;
using namespace Eigen;

// plotter constructor
Plotter::Plotter() {
  // initialize Matlab egine
  _initMatlab();
}

// plotter constructor
Plotter::Plotter(bool toSave, string resultsDir) :
  _toSave(toSave), _resultsDir(resultsDir)
{
  // initialize Matlab egine
  _initMatlab();
  // initialize results directory if asked to save the plots
  if (_toSave) _initResultsDir();
}

// plotter destructor
Plotter::~Plotter() {
//  engClose(_matlab);
}

// initialize Matlab engine
void Plotter::_initMatlab() {
  // open matlab engine
  _matlab = engOpen(NULL);
  // init output buffer
  engOutputBuffer(_matlab, NULL, 0);
}

// initialize results directory
void Plotter::_initResultsDir() {
  mkdir(_resultsDir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
}

// plot the residual map, done with Matlab
void Plotter::plotResidueMap(const MatrixXd& residue,
                             const VectorXd& theta, const VectorXd& phi,
                             double optTheta, double optPhi) {

  // initialize mx pointers to use inside Matlab
  mxArray *mResidue, *mTheta, *mPhi, *mOptTheta, *mOptPhi;

  // theta and phi m-arrays  for residual maps
  mTheta = mxCreateDoubleMatrix(1, theta.size(), mxREAL);
  mPhi = mxCreateDoubleMatrix(1, phi.size(), mxREAL);

  // optimal angle for the residual maps (by GSR)
  mOptTheta = mxCreateDoubleScalar(optTheta);
  mOptPhi = mxCreateDoubleScalar(optPhi);

  // copy data into the pointer destination for the theta and phi m-arrays
  memcpy((double*)mxGetPr(mTheta), (double*)theta.data(),
         theta.size()*sizeof(double));
  memcpy((double*)mxGetPr(mPhi), (double*)phi.data(),
         phi.size()*sizeof(double));

  // put variables inside Matlab
  engPutVariable(_matlab, "theta", mTheta);
  engPutVariable(_matlab, "phi", mPhi);
  engPutVariable(_matlab, "optTheta", mOptTheta);
  engPutVariable(_matlab, "optPhi", mOptPhi);

  // m-array for the residuals
  mResidue = mxCreateDoubleMatrix(residue.rows(), residue.cols(), mxREAL);

  // copy resiuals to Matlab
  memcpy((double*)mxGetPr(mResidue), (double*)residue.data(),
         residue.rows()*residue.cols()*sizeof(double));

  // ... and put it inside Matlab
  engPutVariable(_matlab, "R", mResidue);

  // Matlab command to plot the residual map
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

  // run the command
  engEvalString(_matlab, matlabCommand);
}

// plot the magnetic field map, done with Matlab
void Plotter::plotMagneticMap(const MatrixXd& Axy, const MatrixXd& Bz,
                              const VectorXd& X, const VectorXd& Y,
                              const double Ab) {

  // initialize the pointers to the m-variables
  mxArray *mAxy, *mBz, *mX, *mY, *mAb;

  // initialize boundary vector potential
  mAb = mxCreateDoubleScalar(Ab);

  // initialize m-arrays for the X and Y coordinates in the magnetic field map
  mX = mxCreateDoubleMatrix(1, X.size(), mxREAL);
  mY = mxCreateDoubleMatrix(1, Y.size(), mxREAL);

  // copy the coordinates to the m-arrays
  memcpy((double*)mxGetPr(mX), (double*)X.data(), X.size()*sizeof(double));
  memcpy((double*)mxGetPr(mY), (double*)Y.data(), Y.size()*sizeof(double));

  // put the coordinates into Matlab
  engPutVariable(_matlab, "Ab", mAb);
  engPutVariable(_matlab, "X", mX);
  engPutVariable(_matlab, "Y", mY);

  // create the m-matrix for the vector potential
  mAxy = mxCreateDoubleMatrix(Axy.rows(), Axy.cols(), mxREAL);

  // copy vector potential data into the m-matrix
  memcpy((double*)mxGetPr(mAxy), (double*)Axy.data(),
         Axy.rows()*Axy.cols()*sizeof(double));

  // put it into Matlab
  engPutVariable(_matlab, "Axy", mAxy);

  // create m-matrix for the magnetic field data
  mBz = mxCreateDoubleMatrix(Bz.rows(), Bz.cols(), mxREAL);

  // copy the magnetic field data into the m-matrix
  memcpy((double*)mxGetPr(mBz), (double*)Bz.data(),
         Bz.rows()*Bz.cols()*sizeof(double));

  // ... and put it inside Matlab
  engPutVariable(_matlab, "Bz", mBz);

  // Matlab command for plotting the magnetic field map
  char matlabCommand[] =
  "figure\n"
  "contour(X, Y, Bz, linspace(min(Bz(:)), max(Bz(:)), 100), ...\n"
  "  'Fill', 'on', 'LineStyle', 'none');\n"
  "colorbar\n"
  "hold on\n"
  "contour(X, Y, Axy, 50, 'LineColor', 'black');\n"
  "contour(X, Y, Axy, [Ab Ab], 'LineColor', 'white', 'LineWidth', 3);\n"
  "caxis([min(Bz(:)) max(Bz(:))]);\n"
  "set(gca, 'Color', [0 0 .5]);\n"
  "set(gcf, 'Color', 'white');\n"
  "set(gcf, 'InvertHardCopy', 'off');\n"
  "xlabel 'X_{MC} (AU)'\n"
  "ylabel 'Y_{MC} (AU)'\n";

  // run the final command
  engEvalString(_matlab, matlabCommand);
}

// plot Pt(A) for GSR
void Plotter::plotGsrAPt(const Curve& APtInCurve,
                const Curve& APtOutCurve,
                const Curve& APtFitCurve) {
  Gnuplot APt("gnuplot -persist");
  if (_toSave) {
    APt << "set term post enhanced color eps\n"
        << "set output '" << _resultsDir << "/gsr_Pt_A.eps'\n";
  }
  APt << "p '-' w p t 'Pt(A) in', "
        << "'-' w p t 'Pt(A) out', "
        << "'-' w l t 'Pt(A) fit'\n";
  APt.send(APtInCurve).send(APtOutCurve).send(APtFitCurve);
  if (_toSave) {
    APt << "set term png\n"
        << "set output '" << _resultsDir << "/gsr_Pt_A.png'\n";
    APt << "replot";
  }
}


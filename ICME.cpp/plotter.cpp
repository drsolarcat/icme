
// project headers
#include "plotter.h"
#include "gnuplot.h"
// library headers
#include <Python.h>
#include <numpy/arrayobject.h>
// standard headers
#include <iostream>
#include <string>
#include <sys/stat.h>

using namespace std;
using namespace Eigen;

// plotter constructor
Plotter::Plotter() {
  // initialize Matlab egine
//  _initMatlab();
  // initialize Python
  _initPython();
}

// plotter constructor
Plotter::Plotter(bool toSave, string resultsDir) :
  _toSave(toSave), _resultsDir(resultsDir)
{
  // initialize results directory if asked to save the plots
  if (_toSave) _initResultsDir();
  // initialize Matlab egine
//  _initMatlab();
  // initialize Python
  _initPython();
}

// plotter destructor
Plotter::~Plotter() {
  PyObject_CallObject(
    PyDict_GetItemString(_python_dictionary, "showPlots"), NULL);
//  engClose(_matlab);
}

// initialize Matlab engine
void Plotter::_initMatlab() {
  // open matlab engine
  _matlab = engOpen(NULL);
  // init output buffer
  engOutputBuffer(_matlab, NULL, 0);
}

int Plotter::_initPython() {
  // launch the python interpreter
  Py_Initialize();
  // this macro is defined be NumPy and must be included
  import_array1(-1);
  // update module search path
  PyObject *sys_module = PyImport_ImportModule("sys");
  PyObject *sys_dict = PyModule_GetDict(sys_module);
  PyObject *sys_path = PyMapping_GetItemString(sys_dict, "path");
  PyObject *add_value = PyString_FromString("./");
  PyList_Append(sys_path, add_value);
  Py_DECREF(add_value);
  Py_DECREF(sys_path);
  Py_DECREF(sys_module);
  _python_module = PyImport_ImportModule("plotter");
  _python_dictionary = PyModule_GetDict(_python_module);

  PyObject *pArgs;
  pArgs = PyTuple_New(2);
  PyTuple_SetItem(pArgs, 0, (_toSave ? Py_True : Py_False));
  PyTuple_SetItem(pArgs, 1, PyString_FromString(_resultsDir.c_str()));
  PyObject_CallObject(
    PyDict_GetItemString(_python_dictionary, "initPlotter"), pArgs);
  Py_XDECREF(pArgs);
}

// initialize results directory
void Plotter::_initResultsDir() {
  char strDirEps[256];
  char strDirPng[256];
  strcpy(strDirEps, _resultsDir.c_str());
  strcat(strDirEps, "/eps");
  strcpy(strDirPng, _resultsDir.c_str());
  strcat(strDirPng, "/png");
  mkdir(_resultsDir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
  mkdir(strDirEps, S_IRWXU|S_IRWXG|S_IRWXO);
  mkdir(strDirPng, S_IRWXU|S_IRWXG|S_IRWXO);
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

//  // initialize the pointers to the m-variables
//  mxArray *mAxy, *mBz, *mX, *mY, *mAb;

//  // initialize boundary vector potential
//  mAb = mxCreateDoubleScalar(Ab);

//  // initialize m-arrays for the X and Y coordinates in the magnetic field map
//  mX = mxCreateDoubleMatrix(1, X.size(), mxREAL);
//  mY = mxCreateDoubleMatrix(1, Y.size(), mxREAL);

//  // copy the coordinates to the m-arrays
//  memcpy((double*)mxGetPr(mX), (double*)X.data(), X.size()*sizeof(double));
//  memcpy((double*)mxGetPr(mY), (double*)Y.data(), Y.size()*sizeof(double));

//  // put the coordinates into Matlab
//  engPutVariable(_matlab, "Ab", mAb);
//  engPutVariable(_matlab, "X", mX);
//  engPutVariable(_matlab, "Y", mY);

//  // create the m-matrix for the vector potential
//  mAxy = mxCreateDoubleMatrix(Axy.rows(), Axy.cols(), mxREAL);

//  // copy vector potential data into the m-matrix
//  memcpy((double*)mxGetPr(mAxy), (double*)Axy.data(),
//         Axy.rows()*Axy.cols()*sizeof(double));

//  // put it into Matlab
//  engPutVariable(_matlab, "Axy", mAxy);

//  // create m-matrix for the magnetic field data
//  mBz = mxCreateDoubleMatrix(Bz.rows(), Bz.cols(), mxREAL);

//  // copy the magnetic field data into the m-matrix
//  memcpy((double*)mxGetPr(mBz), (double*)Bz.data(),
//         Bz.rows()*Bz.cols()*sizeof(double));

//  // ... and put it inside Matlab
//  engPutVariable(_matlab, "Bz", mBz);

//  // Matlab command for plotting the magnetic field map
//  char matlabCommand[] =
//  "figure\n"
//  "contour(X, Y, Bz, linspace(min(Bz(:)), max(Bz(:)), 100), ...\n"
//  "  'Fill', 'on', 'LineStyle', 'none');\n"
//  "colorbar\n"
//  "hold on\n"
//  "contour(X, Y, Axy, 50, 'LineColor', 'black');\n"
//  "contour(X, Y, Axy, [Ab Ab], 'LineColor', 'white', 'LineWidth', 3);\n"
//  "caxis([min(Bz(:)) max(Bz(:))]);\n"
//  "set(gca, 'Color', [0 0 .5]);\n"
//  "set(gcf, 'Color', 'white');\n"
//  "set(gcf, 'InvertHardCopy', 'off');\n"
//  "xlabel 'X_{MC} (AU)'\n"
//  "ylabel 'Y_{MC} (AU)'\n";

//  // run the final command
//  engEvalString(_matlab, matlabCommand);

  PyObject *pArgs, *func;

  npy_intp pXDim[] = {X.size()};
  npy_intp pYDim[] = {Y.size()};
  npy_intp pAxyDim[] = {Axy.rows()*Axy.cols()};
  npy_intp pBzDim[] = {Bz.rows()*Bz.cols()};

  pArgs = PyTuple_New(5);

  PyTuple_SetItem(pArgs, 0,
    PyArray_SimpleNewFromData(1, pXDim, PyArray_DOUBLE,
      const_cast<double*>(X.data())));

  PyTuple_SetItem(pArgs, 1,
    PyArray_SimpleNewFromData(1, pYDim, PyArray_DOUBLE,
      const_cast<double*>(Y.data())));

  PyTuple_SetItem(pArgs, 2,
    PyArray_SimpleNewFromData(1, pAxyDim, PyArray_DOUBLE,
      const_cast<double*>(Axy.data())));

  PyTuple_SetItem(pArgs, 3,
    PyArray_SimpleNewFromData(1, pBzDim, PyArray_DOUBLE,
      const_cast<double*>(Bz.data())));

  PyTuple_SetItem(pArgs, 4, PyFloat_FromDouble(Ab));

  func = PyDict_GetItemString(_python_dictionary, "plotBzMap");
  PyObject_CallObject(func, pArgs);

  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot Pt(A) for GSR
void Plotter::plotGsrAPt(const Curve& APtInCurve,
                         const Curve& APtOutCurve,
                         const Curve& APtFitCurve) {

  PyObject *pAPtIn, *pAPtOut, *pAPtFit, *pArgs, *func;

  pAPtIn = PyDict_New();
  npy_intp pAPtInDim[] = {APtInCurve.size()};
  PyDict_SetItemString(pAPtIn, "x",
    PyArray_SimpleNewFromData(1, pAPtInDim, PyArray_DOUBLE,
      const_cast<double*>(APtInCurve.cols().x.data())));
  PyDict_SetItemString(pAPtIn, "y",
    PyArray_SimpleNewFromData(1, pAPtInDim, PyArray_DOUBLE,
      const_cast<double*>(APtInCurve.cols().y.data())));

  pAPtOut = PyDict_New();
  npy_intp pAPtOutDim[] = {APtOutCurve.size()};
  PyDict_SetItemString(pAPtOut, "x",
    PyArray_SimpleNewFromData(1, pAPtOutDim, PyArray_DOUBLE,
      const_cast<double*>(APtOutCurve.cols().x.data())));
  PyDict_SetItemString(pAPtOut, "y",
    PyArray_SimpleNewFromData(1, pAPtOutDim, PyArray_DOUBLE,
      const_cast<double*>(APtOutCurve.cols().y.data())));

  pAPtFit = PyDict_New();
  npy_intp pAPtFitDim[] = {APtFitCurve.size()};
  PyDict_SetItemString(pAPtFit, "x",
    PyArray_SimpleNewFromData(1, pAPtFitDim, PyArray_DOUBLE,
      const_cast<double*>(APtFitCurve.cols().x.data())));
  PyDict_SetItemString(pAPtFit, "y",
    PyArray_SimpleNewFromData(1, pAPtFitDim, PyArray_DOUBLE,
      const_cast<double*>(APtFitCurve.cols().y.data())));

  pArgs = PyTuple_New(3);
  PyTuple_SetItem(pArgs, 0, pAPtIn);
  PyTuple_SetItem(pArgs, 1, pAPtOut);
  PyTuple_SetItem(pArgs, 2, pAPtFit);
  func = PyDict_GetItemString(_python_dictionary, "plotGsrAPt");
  PyObject_CallObject(func, pArgs);

  Py_XDECREF(pAPtIn);
  Py_XDECREF(pAPtOut);
  Py_XDECREF(pAPtFit);
  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot dPt/dA(A) for GSR
void Plotter::plotGsrAdPt(const Curve& AdPtFitCurve) {
  Gnuplot AdPt("gnuplot -persist");
  if (_toSave) {
    AdPt << "set term post enhanced color eps\n"
        << "set output '" << _resultsDir << "/gsr_A_dPt.eps'\n";
  }
  AdPt << "p '-' w l t 'dPt/dA fit'\n";
  AdPt.send(AdPtFitCurve);
}

// plot Bz(A) for GSR
void Plotter::plotGsrABz(const Curve& ABzCurve, const Curve& ABzFitCurve) {
  Gnuplot ABz("gnuplot -persist");
  if (_toSave) {
    ABz << "set term post enhanced color eps\n"
        << "set output '" << _resultsDir << "/gsr_A_Bz.eps'\n";
  }
  ABz << "p '-' w p t 'Bz(A)', "
      << "'-' w l t 'Bz(A) fit'\n";
  ABz.send(ABzCurve).send(ABzFitCurve);
}


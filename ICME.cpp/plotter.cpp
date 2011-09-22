
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
  // initialize Python
  _initPython();
}

// plotter constructor
Plotter::Plotter(bool toSave, string resultsDir) :
  _toSave(toSave), _resultsDir(resultsDir)
{
  // initialize results directory if asked to save the plots
  if (_toSave) _initResultsDir();
  // initialize Python
  _initPython();
}

// plotter destructor
Plotter::~Plotter() {
  PyObject_CallObject(
    PyDict_GetItemString(_python_dictionary, "showPlots"), NULL);
//  engClose(_matlab);
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

  PyObject *pArgs, *func;

  npy_intp pThetaDim[] = {theta.size()};
  npy_intp pPhiDim[] = {phi.size()};
  npy_intp pResidueDim[] = {residue.rows()*residue.cols()};

  pArgs = PyTuple_New(3);

  PyTuple_SetItem(pArgs, 0,
    PyArray_SimpleNewFromData(1, pThetaDim, PyArray_DOUBLE,
      const_cast<double*>(theta.data())));

  PyTuple_SetItem(pArgs, 1,
    PyArray_SimpleNewFromData(1, pPhiDim, PyArray_DOUBLE,
      const_cast<double*>(phi.data())));

  PyTuple_SetItem(pArgs, 2,
    PyArray_SimpleNewFromData(1, pResidueDim, PyArray_DOUBLE,
      const_cast<double*>(residue.data())));

  func = PyDict_GetItemString(_python_dictionary, "plotGsrResidue");
  PyObject_CallObject(func, pArgs);

  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot the magnetic field map, done with Matlab
void Plotter::plotMagneticMap(const MatrixXd& Axy, const MatrixXd& Bz,
                              const VectorXd& X, const VectorXd& Y,
                              const double Ab) {

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



// project headers
#include "plotter.h"
// library headers
#include <Python.h>
#include <numpy/arrayobject.h>
// standard headers
#include <string>
#include <sys/stat.h>

using namespace std;
using namespace Eigen;

// plotter constructor
Plotter::Plotter()
{
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
Plotter::~Plotter()
{
  // show all pending matplotlib plots
  PyObject_CallObject(
    PyDict_GetItemString(_python_dictionary, "showPlots"), NULL);
}

int Plotter::_initPython()
{
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
  // remove temporary references
  Py_DECREF(add_value);
  Py_DECREF(sys_path);
  Py_DECREF(sys_module);
  // assign python and dictionary objects
  _python_module = PyImport_ImportModule("plotter");
  _python_dictionary = PyModule_GetDict(_python_module);

  // initializing the matplotlib plotter
  PyObject *pArgs; // arguments for a python function
  pArgs = PyTuple_New(2); // initialize the arguments tuple
  // whether to save the plots to files
  PyTuple_SetItem(pArgs, 0, (_toSave ? Py_True : Py_False));
  // directory where to save the plots
  PyTuple_SetItem(pArgs, 1, PyString_FromString(_resultsDir.c_str()));
  // call the matplotlib initialization script
  PyObject_CallObject(
    PyDict_GetItemString(_python_dictionary, "initPlotter"), pArgs);
  // remove reference to arguments tuple
  Py_XDECREF(pArgs);
}

// initialize results directory
void Plotter::_initResultsDir()
{
  char strDirEps[256]; // cstring for eps plots directory
  char strDirPng[256]; // cstring for png plots directory
  // construct eps directory string
  strcpy(strDirEps, _resultsDir.c_str());
  strcat(strDirEps, "/eps");
  // construct png directory string
  strcpy(strDirPng, _resultsDir.c_str());
  strcat(strDirPng, "/png");
  // create plots directory
  mkdir(_resultsDir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
  // create eps plots subdirectory
  mkdir(strDirEps, S_IRWXU|S_IRWXG|S_IRWXO);
  // create png plots subdirectory
  mkdir(strDirPng, S_IRWXU|S_IRWXG|S_IRWXO);
}

// plot the residual map for GSR
void Plotter::plotGsrResidueMap(const MatrixXd& residue,
                                const VectorXd& theta, const VectorXd& phi,
                                double optTheta, double optPhi,
                                double mvabTheta, double mvabPhi,
                                double mvubTheta, double mvubPhi)
{
  PyObject *pArgs, *func; // pointers to arguments and function objects

  pArgs = PyTuple_New(mvabTheta && mvabPhi ? 7 : 5); // initialize the arguments tuple

  // set theta array argument
  npy_intp pThetaDim[] = {theta.size()}; // shape array
  PyTuple_SetItem(pArgs, 0,
    PyArray_SimpleNewFromData(1, pThetaDim, PyArray_DOUBLE,
      const_cast<double*>(theta.data())));

  // set phi array argument
  npy_intp pPhiDim[] = {phi.size()}; // shape array
  PyTuple_SetItem(pArgs, 1,
    PyArray_SimpleNewFromData(1, pPhiDim, PyArray_DOUBLE,
      const_cast<double*>(phi.data())));

  // set residue array argument
  // here we pass the matrix as a vector and will reshape it later in python
  npy_intp pResidueDim[] = {residue.rows()*residue.cols()}; // shape array
  PyTuple_SetItem(pArgs, 2,
    PyArray_SimpleNewFromData(1, pResidueDim, PyArray_DOUBLE,
      const_cast<double*>(residue.data())));

  // set optimal theta argument
  PyTuple_SetItem(pArgs, 3, PyFloat_FromDouble(optTheta));

  // set optimal phi argument
  PyTuple_SetItem(pArgs, 4, PyFloat_FromDouble(optPhi));

  if (mvabTheta && mvabPhi) {
    // set MVAB theta argument
    PyTuple_SetItem(pArgs, 5, PyFloat_FromDouble(mvabTheta));
    // set MVAB phi argument
    PyTuple_SetItem(pArgs, 6, PyFloat_FromDouble(mvabPhi));
  }

  // initialize and call the python function
  func = PyDict_GetItemString(_python_dictionary, "plotGsrResidue");
  PyObject_CallObject(func, pArgs);
}

// plot the magnetic field map for GSR
void Plotter::plotGsrMagneticMap(const MatrixXd& Axy, const MatrixXd& Bz,
                                 const VectorXd& X, const VectorXd& Y,
                                 const double Ab)
{
  PyObject *pArgs, *func; // pointers to the arguments and function objects

  // initialize the shape arrays
  npy_intp pXDim[] = {X.size()};
  npy_intp pYDim[] = {Y.size()};
  // will pass matrix as vector and reshape it later
  npy_intp pAxyDim[] = {Axy.rows()*Axy.cols()};
  // will pass matrix as vector and reshape it later
  npy_intp pBzDim[] = {Bz.rows()*Bz.cols()};

  pArgs = PyTuple_New(5); // initialize the arguments tuple

  // set X
  PyTuple_SetItem(pArgs, 0,
    PyArray_SimpleNewFromData(1, pXDim, PyArray_DOUBLE,
      const_cast<double*>(X.data())));

  // set Y
  PyTuple_SetItem(pArgs, 1,
    PyArray_SimpleNewFromData(1, pYDim, PyArray_DOUBLE,
      const_cast<double*>(Y.data())));

  // set Axy
  PyTuple_SetItem(pArgs, 2,
    PyArray_SimpleNewFromData(1, pAxyDim, PyArray_DOUBLE,
      const_cast<double*>(Axy.data())));

  // set Bz
  PyTuple_SetItem(pArgs, 3,
    PyArray_SimpleNewFromData(1, pBzDim, PyArray_DOUBLE,
      const_cast<double*>(Bz.data())));

  // set boundary vector potential Ab
  PyTuple_SetItem(pArgs, 4, PyFloat_FromDouble(Ab));

  // initialize and call the python function
  func = PyDict_GetItemString(_python_dictionary, "plotGsrBzMap");
  PyObject_CallObject(func, pArgs);

  // delete the references
  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot Pt(A) for GSR
void Plotter::plotGsrAPt(const Curve& APtInCurve,
                         const Curve& APtOutCurve,
                         const Curve& APtFitCurve)
{
  // pointers to the python objects
  PyObject *pAPtIn, *pAPtOut, *pAPtFit, *pArgs, *func;

  // set inward Pt(A)
  pAPtIn = PyDict_New(); // tuple for [x,y] coordinates
  npy_intp pAPtInDim[] = {APtInCurve.size()}; // shape array
  // set x
  PyDict_SetItemString(pAPtIn, "x",
    PyArray_SimpleNewFromData(1, pAPtInDim, PyArray_DOUBLE,
      const_cast<double*>(APtInCurve.cols().x.data())));
  // set y
  PyDict_SetItemString(pAPtIn, "y",
    PyArray_SimpleNewFromData(1, pAPtInDim, PyArray_DOUBLE,
      const_cast<double*>(APtInCurve.cols().y.data())));

  // set outward Pt(A)
  pAPtOut = PyDict_New();
  npy_intp pAPtOutDim[] = {APtOutCurve.size()};
  PyDict_SetItemString(pAPtOut, "x",
    PyArray_SimpleNewFromData(1, pAPtOutDim, PyArray_DOUBLE,
      const_cast<double*>(APtOutCurve.cols().x.data())));
  PyDict_SetItemString(pAPtOut, "y",
    PyArray_SimpleNewFromData(1, pAPtOutDim, PyArray_DOUBLE,
      const_cast<double*>(APtOutCurve.cols().y.data())));

  // set fitted Pt(A)
  pAPtFit = PyDict_New();
  npy_intp pAPtFitDim[] = {APtFitCurve.size()};
  PyDict_SetItemString(pAPtFit, "x",
    PyArray_SimpleNewFromData(1, pAPtFitDim, PyArray_DOUBLE,
      const_cast<double*>(APtFitCurve.cols().x.data())));
  PyDict_SetItemString(pAPtFit, "y",
    PyArray_SimpleNewFromData(1, pAPtFitDim, PyArray_DOUBLE,
      const_cast<double*>(APtFitCurve.cols().y.data())));

  pArgs = PyTuple_New(3); // tuple for arguments
  // fill in the tuple
  PyTuple_SetItem(pArgs, 0, pAPtIn);
  PyTuple_SetItem(pArgs, 1, pAPtOut);
  PyTuple_SetItem(pArgs, 2, pAPtFit);

  // initialize and call the python function
  func = PyDict_GetItemString(_python_dictionary, "plotGsrAPt");
  PyObject_CallObject(func, pArgs);

  // delete the references
  Py_XDECREF(pAPtIn);
  Py_XDECREF(pAPtOut);
  Py_XDECREF(pAPtFit);
  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot dPt/dA(A) for GSR
void Plotter::plotGsrAdPt(const Curve& AdPtFitCurve)
{
  PyObject *pAdPt, *pArgs, *func; // pointers to python objects

  // set dPt/dA
  pAdPt = PyDict_New();
  npy_intp pAdPtDim[] = {AdPtFitCurve.size()}; // shape
  // x
  PyDict_SetItemString(pAdPt, "x",
    PyArray_SimpleNewFromData(1, pAdPtDim, PyArray_DOUBLE,
      const_cast<double*>(AdPtFitCurve.cols().x.data())));
  // y
  PyDict_SetItemString(pAdPt, "y",
    PyArray_SimpleNewFromData(1, pAdPtDim, PyArray_DOUBLE,
      const_cast<double*>(AdPtFitCurve.cols().y.data())));

  // initialize and fill in the arguments tuple
  pArgs = PyTuple_New(1);
  PyTuple_SetItem(pArgs, 0, pAdPt);

  // initialize and call the python function
  func = PyDict_GetItemString(_python_dictionary, "plotGsrAdPt");
  PyObject_CallObject(func, pArgs);

  // delete the references
  Py_XDECREF(pAdPt);
  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot Bz(A) for GSR
void Plotter::plotGsrABz(const Curve& ABzCurve, const Curve& ABzFitCurve)
{
  PyObject *pABz, *pABzFit, *pArgs, *func; // pointers to the python objects

  // set Bz(A)
  pABz = PyDict_New();
  npy_intp pABzDim[] = {ABzCurve.size()}; // shape
  PyDict_SetItemString(pABz, "x",
    PyArray_SimpleNewFromData(1, pABzDim, PyArray_DOUBLE,
      const_cast<double*>(ABzCurve.cols().x.data())));
  PyDict_SetItemString(pABz, "y",
    PyArray_SimpleNewFromData(1, pABzDim, PyArray_DOUBLE,
      const_cast<double*>(ABzCurve.cols().y.data())));

  // set fitted Bz(A)
  pABzFit = PyDict_New();
  npy_intp pABzFitDim[] = {ABzFitCurve.size()};
  PyDict_SetItemString(pABzFit, "x",
    PyArray_SimpleNewFromData(1, pABzFitDim, PyArray_DOUBLE,
      const_cast<double*>(ABzFitCurve.cols().x.data())));
  PyDict_SetItemString(pABzFit, "y",
    PyArray_SimpleNewFromData(1, pABzFitDim, PyArray_DOUBLE,
      const_cast<double*>(ABzFitCurve.cols().y.data())));

  // initialize and fill in the arguments tuple
  pArgs = PyTuple_New(2);
  PyTuple_SetItem(pArgs, 0, pABz);
  PyTuple_SetItem(pArgs, 1, pABzFit);

  // initialize and call the python function
  func = PyDict_GetItemString(_python_dictionary, "plotGsrABz");
  PyObject_CallObject(func, pArgs);

  // delete the references
  Py_XDECREF(pABz);
  Py_XDECREF(pABzFit);
  Py_XDECREF(pArgs);
  Py_XDECREF(func);
}

// plot B rotation for MVA
void Plotter::plotMvaBrot(const Eigen::VectorXd& Bx, const Eigen::VectorXd& By)
{
  PyObject *pArgs, *func; // pointers to the python object

  // initialize the shape arrays
  npy_intp pBxDim[] = {Bx.size()};
  npy_intp pByDim[] = {By.size()};

  pArgs = PyTuple_New(2); // initialize the arguments tuple

  // set Bx
  PyTuple_SetItem(pArgs, 0,
    PyArray_SimpleNewFromData(1, pBxDim, PyArray_DOUBLE,
      const_cast<double*>(Bx.data())));

  // set By
  PyTuple_SetItem(pArgs, 1,
    PyArray_SimpleNewFromData(1, pByDim, PyArray_DOUBLE,
      const_cast<double*>(By.data())));

  // initialize and call the python function
  func = PyDict_GetItemString(_python_dictionary, "plotMvaBrot");
  PyObject_CallObject(func, pArgs);
}


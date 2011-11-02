
// project headers
#include "curve.h"
#include "event.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <Python.h>
// standard headers
#include <string>

// a class for drawing different kinds of plots
class Plotter {
  PyObject* _python_module; // python module
  PyObject* _python_dictionary; // python dictionary
  bool _toSave; // whether to save the plot or not
  std::string _resultsDir; // directory for resulting plots to be saved into
  public:
    // constructors
    Plotter();
    Plotter(bool, std::string);
    ~Plotter(); // destructor
    // accessors
    const bool toSave() const {return _toSave;}
    const std::string resultsDir() const {return _resultsDir;}
    // plot residue map
    void plotGsrResidueMap(const Eigen::MatrixXd&,
                           const Eigen::VectorXd&, const Eigen::VectorXd&,
                           double, double,
                           double mvabTheta=NULL, double mvabPhi=NULL,
                           double mvubTheta=NULL, double mvubPhi=NULL);
    // plot magnetic field map
    void plotGsrMagneticMap(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                            const Eigen::VectorXd&, const Eigen::VectorXd&,
                            const double, const double);
    // plot Pt(A) plot for GSR
    void plotGsrAPt(const Curve&, const Curve&, const Curve&);
    // plot dPt/dA(A) plot for GSR
    void plotGsrAdPt(const Curve&);
    // plot Bz(A) plot for GSR
    void plotGsrABz(const Curve&, const Curve&);
    // plot B rotation for MVA
    void plotMvaBrot(const Eigen::VectorXd&, const Eigen::VectorXd&);
    // plot all the in-situ data
    void plotData(const Event&);
  protected:
    // initialize Python
    int _initPython();
    // initialize results directory
    void _initResultsDir();
};


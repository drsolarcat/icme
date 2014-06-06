
// project headers
#include "curve.h"
#include "event.h"
#include "axes.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <python2.7/Python.h>
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
                           double mvabTheta=0, double mvabPhi=0,
                           double mvubTheta=0, double mvubPhi=0);
    // plot magnetic field map
    void plotGsrMagneticMap(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                            const Eigen::VectorXd&, const Eigen::VectorXd&,
                            const double, const double,
                            const Eigen::VectorXd&, const Eigen::VectorXd&,
                            const Axes&);
    // plot Pt(A) plot for GSR
    void plotGsrAPt(const Curve&, const Curve&, const Curve&);
    // plot Pt(A) plot for GSR for all the analyzed time period
    void plotGsrAPtFull(const Curve&);
    // plot dPt/dA(A) plot for GSR
    void plotGsrAdPt(const Curve&);
    // plot Bz(A) plot for GSR
    void plotGsrABz(const Curve&, const Curve&);
    // plot B rotation for MVA
    void plotMvaBrot(const Eigen::VectorXd&, const Eigen::VectorXd&);
    // plot all the in-situ data
    void plotData(const Event&);
    // plot B rotation
    void plotBrot(const Event&);
    // plot simple 1D data
    void plotData1D(const Eigen::VectorXd&);
  protected:
    // initialize Python
    int _initPython();
    // initialize results directory
    void _initResultsDir();
};


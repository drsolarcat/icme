
// project headers
#include "curve.h"
// library headers
#include "engine.h"
#include <eigen3/Eigen/Dense>
#include <Python.h>
// standard headers
#include <string>

// a class for drawing different kinds of plots
class Plotter {
  Engine* _matlab; // matlab engine object
  PyObject* _python_module; // python module
  PyObject* _python_dictionary; // python dictionary
  public:
    // constructors
    Plotter();
    Plotter(bool, std::string);
    ~Plotter(); // destructor
    bool _toSave; // whether to save the plot or not
    std::string _resultsDir; // directory for resulting plots to be saved into
    // accessors
    const bool toSave() const {return _toSave;}
    const std::string resultsDir() const {return _resultsDir;}
    // plot residue map
    void plotResidueMap(const Eigen::MatrixXd&,
                        const Eigen::VectorXd&, const Eigen::VectorXd&,
                        double, double);
    // plot magnetic field map
    void plotMagneticMap(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                         const Eigen::VectorXd&, const Eigen::VectorXd&,
                         const double);
    // plot Pt(A) plot for GSR
    void plotGsrAPt(const Curve&, const Curve&, const Curve&);
    // plot dPt/dA(A) plot for GSR
    void plotGsrAdPt(const Curve&);
    // plot Bz(A) plot for GSR
    void plotGsrABz(const Curve&, const Curve&);
  protected:
    // initialize Matlab engine
    void _initMatlab();
    // initialize Python
    int _initPython();
    // initialize results directory
    void _initResultsDir();
};


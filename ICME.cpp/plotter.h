
// library headers
#include "engine.h"
#include <eigen3/Eigen/Dense>

// a class for drawing different kinds of plots
class Plotter {
  Engine* _matlab;
  public:
    Plotter();
    ~Plotter();
    void plotResidueMap(const Eigen::MatrixXd&,
                        const Eigen::VectorXd&, const Eigen::VectorXd&,
                        double, double);
    void plotMagneticMap(const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                         const Eigen::VectorXd&, const Eigen::VectorXd&,
                         const double);
};


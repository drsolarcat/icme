
#ifndef MVA_ANALYZER_H
#define MVA_ANALYZER_H

#include "event.h"

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

// this class is used for MVA analysis
class MvaAnalyzer {
  public:
    void analyze(Event&); // main trigger to carry the analysis,
                          // the results are stored in the Event object
    // ordinary minimum variance analysis
    void analyzeMvab(Event&);
    // minimum variance analysis with unit (normalized) vectors
    void analyzeMvub(Event&);
  private:
    // calculate variance matrix from three vectors of data (x, y, z)
    Matrix3d calculateVarianceMatrix(const VectorXd&,
                                     const VectorXd&,
                                     const VectorXd&);
    MvaResults analyzeMva(const Event&, bool);
};

#endif


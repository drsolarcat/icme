
#ifndef MVA_ANALYZER_H
#define MVA_ANALYZER_H

#include "event.h"

#include <eigen3/Eigen/Dense>

// this class is used for MVA analysis
class MvaAnalyzer {
  public:
    void analyze(Event&); // main trigger to carry the analysis,
                          // the results are stored in the Event object
    // ordinary minimum variance analysis
    void analyzeMvab(Event&);
    void analyzePmvab(Event&); // projected
    // minimum variance analysis with unit (normalized) vectors
    void analyzeMvub(Event&);
    void analyzePmvub(Event&); // projected
  private:
    // calculate variance matrix from three vectors of data (x, y, z)
    Eigen::Matrix3d calculateVarianceMatrix(const Eigen::VectorXd&,
                                            const Eigen::VectorXd&,
                                            const Eigen::VectorXd&);
    // MVA analysis
    MvaResults analyzeMva(const Event&, bool);
    // PMVA analysis
    PmvaResults analyzePmva(const Event&, bool);
};

#endif


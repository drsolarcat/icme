
#ifndef EVENT_H
#define EVENT_H

#include "data.h"
#include "config.h"
#include "axes.h"

#include <eigen3/Eigen/Dense>

#include <vector>

// this structure holds the results of MVA analysis
struct MvaResults {
  Axes axes; //axes
  Eigen::Matrix3d matrix; // variance matrix
  double lambdaMin, lambdaMed, lambdaMax; // eigen values
  double criterion; // variance analysis criterion
};

// this structure holds the results of projected MVA analysis
struct PmvaResults {
  Axes axes; //axes
  Eigen::Matrix3d matrix, // variance matrix = M
                  pmatrix, // projection matrix = P
                  pmpmatrix; // projected variance matrix = PMP
  double lambdaMin, lambdaMed, lambdaMax; // eigen values of PMP
};

// this structure holds results of deHoffmann-Teller analysis
struct DhtResults {
  Eigen::Vector3d Vht; // deHoffmann-Teller frame speed
  double cc; // correlation coefficient of dHT analysis
};

// this structure holds results of a single run of axes search algorithm
struct GsrRun {
  double minTheta, dTheta, maxTheta, // theta limits
         minPhi, dPhi, maxPhi, // phi limits
         optTheta, optPhi; // optimal theta and phi
  Eigen::MatrixXd originalResidue, // original residue
                  combinedResidue, // renormalized residue
                  branchLength; // branch length
};

// this structure holds results of GSR analysis
struct GsrResults {
  std::vector<GsrRun> runs; // runs of axes estimation loop
};

// a class representing one event, it stores all the data and results of
// analysis of one particular event
class Event {
    Data e_dataWide, e_dataNarrow; // data
    ConfigRow e_config; // configuration structure
    // results of MVAB and MVUB analysis
    MvaResults e_mvab, e_mvub;
    // results of projected MVAB and MVUB analysis
    PmvaResults e_pmvab, e_pmvub;
    GsrResults e_gsr; // results of GSR analysis
    DhtResults e_dht; // results of dHT analysis
//    bool e_isMva, e_isDht, e_isGsr;
  public:
    Event(ConfigRow, Data, Data); // construct the initial object with data
                                  // and configuration
    const ConfigRow& config() const {return e_config;} // config getter
    const Data& dataWide() const {return e_dataWide;} // dataWide getter
    const Data& dataNarrow() const {return e_dataNarrow;} // dataNarrow getter
    // MVA results setters and getters
    Event& mvab(MvaResults mvab) {e_mvab = mvab; return *this;}
    const MvaResults& mvab() const {return e_mvab;}
    Event& mvub(MvaResults mvub) {e_mvub = mvub; return *this;}
    const MvaResults& mvub() const {return e_mvub;}
    Event& pmvab(PmvaResults pmvab) {e_pmvab = pmvab; return *this;}
    const PmvaResults& pmvab() const {return e_pmvab;}
    Event& pmvub(PmvaResults pmvub) {e_pmvub = pmvub; return *this;}
    const PmvaResults& pmvub() const {return e_pmvub;}
    // dHT setter and getter
    Event& dht(DhtResults dht) {e_dht = dht; return *this;}
    const DhtResults& dht() const {return e_dht;}
    // GSR setter and getter
    Event& gsr(GsrResults gsr) {e_gsr = gsr; return *this;}
    const GsrResults& gsr() const {return e_gsr;}
};

#endif


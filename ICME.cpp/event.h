
#ifndef EVENT_H
#define EVENT_H

//project headers
#include "data.h"
#include "config.h"
#include "axes.h"
#include "branched_curve.h"
#include "curve.h"
// library headers
#include <eigen3/Eigen/Dense>
// standard headers
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

// this structure holds results of GSR analysis
struct GsrResults {
  double minTheta, dTheta, maxTheta, // theta limits
         minPhi, dPhi, maxPhi, // phi limits
         optTheta, optPhi; // optimal theta and phi
  Eigen::MatrixXd originalResidue, // original residue
                  combinedResidue, // renormalized residue
                  branchLength; // branch length
  Axes axes; // MC axes
  BranchedCurve curve; // Pt(A) curve
  double dx, dy; // steps used for reconstruction
  double Ab, Ac; // boundary and central values of the vector potential
  int Nx, Ny; // number of steps
  Eigen::VectorXd X, Y; // reconstruction coordinates
  Eigen::MatrixXd Axy, Bz; // potential and magnetic field maps
  Curve APtInCurve, APtOutCurve, // inward and outward branches of Pt(A)
        APtFitCurve, AdPtFitCurve, // Pt(A) and dPt/dA(A) fits
        ABzCurve, ABzFitCurve; // Bz(A) curve and its fit
};

// a class representing one event, it stores all the data and results of
// analysis of one particular event
class Event {
    Data _dataWide, _dataNarrow; // data
    ConfigRow _config; // configuration structure
    // results of MVAB and MVUB analysis
    MvaResults _mvab, _mvub;
    // results of projected MVAB and MVUB analysis
    PmvaResults _pmvab, _pmvub;
    GsrResults _gsr; // results of GSR analysis
    DhtResults _dht; // results of dHT analysis
//    bool e_isMva, e_isDht, e_isGsr;
  public:
    Event(ConfigRow, Data, Data); // construct the initial object with data
                                  // and configuration
    const ConfigRow& config() const {return _config;} // config getter
    const Data& dataWide() const {return _dataWide;} // dataWide getter
    const Data& dataNarrow() const {return _dataNarrow;} // dataNarrow getter
    // MVA results setters and getters
    Event& mvab(MvaResults mvab) {_mvab = mvab; return *this;}
    const MvaResults& mvab() const {return _mvab;}
    Event& mvub(MvaResults mvub) {_mvub = mvub; return *this;}
    const MvaResults& mvub() const {return _mvub;}
    Event& pmvab(PmvaResults pmvab) {_pmvab = pmvab; return *this;}
    const PmvaResults& pmvab() const {return _pmvab;}
    Event& pmvub(PmvaResults pmvub) {_pmvub = pmvub; return *this;}
    const PmvaResults& pmvub() const {return _pmvub;}
    // dHT setter and getter
    Event& dht(DhtResults dht) {_dht = dht; return *this;}
    const DhtResults& dht() const {return _dht;}
    // GSR setter and getter
    Event& gsr(GsrResults gsr) {_gsr = gsr; return *this;}
    const GsrResults& gsr() const {return _gsr;}
};

#endif


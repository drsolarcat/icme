
#ifndef EVENT_H
#define EVENT_H

#include <eigen3/Eigen/Dense>

#include "data.h"
#include "config.h"

struct MvaResults {
  Eigen::Vector3d x, y, z;
  Eigen::Matrix3d matrix;
  double lambdaMin, lambdaMed, lambdaMax;
  double criterion;
};

struct DhtResults {
  Eigen::Vector3d Vht;
  double cc;
};

struct GsrResults {

};

// a class representing one event, it stores all the data and results of
// analysis of one particular event
class Event {
    Data e_dataWide, e_dataNarrow; // data
    ConfigRow e_config; // configuration structure
    MvaResults e_mvab, e_mvub; // results of MVAB and MVUB analysis
    GsrResults e_gsr; // results of GSR analysis
    DhtResults e_dht; // results of dHT analysis
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
    // dHT setter and getter
    Event& dht(DhtResults dht) {e_dht = dht; return *this;}
    const DhtResults& dht() const {return e_dht;}
    // GSR setter and getter
    Event& gsr(GsrResults gsr) {e_gsr = gsr; return *this;}
    const GsrResults& gsr() const {return e_gsr;}
};

#endif


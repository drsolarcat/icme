
#ifndef EVENT_H
#define EVENT_H

#include <eigen3/Eigen/Dense>

#include "data.h"
#include "config.h"

using namespace std;
using namespace Eigen;

struct MvaResults {
  Vector3d x, y, z;
  Matrix3d matrix;
  double lambdaMin, lambdaMed, lambdaMax;
  double criterion;
};

// a class representing one event, it stores all the data and results of
// analysis of one particular event
class Event {
    Data e_dataWide, e_dataNarrow; // data
    ConfigRow e_config; // configuration structure
    MvaResults e_mvab, e_mvub;
  public:
    Event(ConfigRow, Data, Data); // construct the initial object with data
                                  // and configuration
    const ConfigRow& config() const {return e_config;}
    const Data& dataWide() const {return e_dataWide;}
    const Data& dataNarrow() const {return e_dataNarrow;}
    Event& mvab(MvaResults mvab) {e_mvab = mvab; return *this;}
    const MvaResults& mvab() const {return e_mvab;}
    Event& mvub(MvaResults mvub) {e_mvub = mvub; return *this;}
    const MvaResults& mvub() const {return e_mvub;}
};

#endif


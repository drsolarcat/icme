
#include "gsr_curve.h"
#include "axes.h"

#include <gsl/gsl_const_mksa.h>

using namespace std;
using namespace Eigen;

// construct the Pt(A) curve
GsrCurve::GsrCurve(Event& event, Axes axes) {
  // step in x direction (sunword)
  const double dx = -event.dht().Vht.dot(axes.x)*240;

  Data data(event.dataNarrow()); // copy data

  data.project(axes); // project data to temporary axes

  // iterate through data
  for (int i=0; i < data.rows().size(); i++) {
    c_vectors.x(i) = (i == 0) ? 0 : 1;
//    if (i == 0) { // the first point of vector potential is 0
//      c_vectors.x(i) = 0;
//    } else { // otherwise perform numerical integration
//      c_vectors.x(i) = 1; // TODO
//    }
    // transverse pressure
    c_vectors.y(i) = data.row(i).Pth+
                     pow(data.row(i).Bz, 2)/2/
                     GSL_CONST_MKSA_VACUUM_PERMEABILITY;
  } // end iteration through the data
}


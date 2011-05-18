
#include "gsr_curve.h"
#include "axes.h"
#include "integrator.h"

#include <gsl/gsl_const_mksa.h>

using namespace std;
using namespace Eigen;

// construct the Pt(A) curve
GsrCurve::GsrCurve(Event& event, Axes axes) {
  // step in x direction (sunword)
  const double dx = -event.dht().Vht.dot(axes.x)*
                     event.config().samplingInterval;

  Integrator integrator;

  Data data(event.dataNarrow()); // copy data

  const int n = data.rows().size();

  data.project(axes); // project data to temporary axes

  _vectors.x = VectorXd::Zero(n);
  _vectors.y = VectorXd::Zero(n);

  // iterate through data
  for (int i=0; i < n; i++) {
    _vectors.x(i) = (i == 0) ? 0 : 1;
    if (i == 0) { // the first point of vector potential is 0
      _vectors.x(i) = 0;
    } else { // otherwise perform numerical integration
      _vectors.x(i) = integrator.NewtonCotes(4,
                                             VectorXd::LinSpaced(i+1, 0, dx*i),
                                             data.cols().By.head(i+1));
    }
    // transverse pressure
    _vectors.y(i) = data.row(i).Pth+
                    pow(data.row(i).Bz, 2)/2/
                    GSL_CONST_MKSA_VACUUM_PERMEABILITY;
  } // end iteration through the data
}


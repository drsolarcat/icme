
// project headers
#include "gsr_analyzer.h"
#include "dht_analyzer.h"
#include "mva_analyzer.h"
#include "gsr_curve.h"
#include "fit_poly.h"
#include "fit_exp.h"
#include "fit_poly_exp.h"
#include "differentiator.h"
#include "my_time.h"
#include "plotter.h"
#include "filter.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_fit.h>
#include <log4cplus/logger.h>
#include <log4cplus/configurator.h>
#include <boost/algorithm/string.hpp>
#include <cxform/cxform.h>
// standard headers
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;
using namespace log4cplus;
using namespace My;

// main trigger to perform GSR analysis of the event
void GsrAnalyzer::analyze(Event& event) {
  DhtAnalyzer dht; // initialize dHT analyzer
  MvaAnalyzer mva; // initialize MVA analyzer

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  LOG4CPLUS_DEBUG(logger, "starting dHT analysis");
  dht.analyze(event); // carry dHT analysis for the event

  LOG4CPLUS_DEBUG(logger, "starting PMVAB analysis to find initial axes");
  mva.analyzePmvab(event); // carry projected MVA anaysis to get initial axes

  // replace initial axes with GSE(RTN) axes
  /*
  PmvaResults pmva;
  Axes pmvaAxes;
  pmvaAxes.x = Vector3d(1, 0, 0);
  pmvaAxes.y = Vector3d(0, 0 ,1);
  pmvaAxes.z = Vector3d(0, -1, 0);
  pmva.axes = pmvaAxes;
  event.pmvab(pmva);
  */

  // make a run of axes searching algorithm, save the results in the gsr
  // structure
  LOG4CPLUS_DEBUG(logger, "searching for optimal axes with 1 degree step");
  GsrResults gsr = loopAxes(event, 0, 1, 90, 0, 1, 360);

  detectAxes(event, gsr);
//  detectAxesByResidue(event, gsr);
//  detectAxesByVariance(event, gsr);

  // flip axes if necessary
  Data dataTmp(event.dataNarrow());
  dataTmp.project(gsr.axes);
  VectorXd xTmp = VectorXd::LinSpaced(dataTmp.rows().size(), 0,
                                      dataTmp.rows().size()-1);
  PolyFit BzTmpFit(dataTmp.rows().size(), xTmp.data(),
                   const_cast<double*>(dataTmp.cols().Bz.data()), 2);
  BzTmpFit.fit();
  if (BzTmpFit.c()[2] > 0) {
    std::cout << "reversed" << std::endl;
    gsr.axes.z = -gsr.axes.z;
    gsr.axes.y = -gsr.axes.y;
  }

  LOG4CPLUS_DEBUG(logger, "magnetic cloud axes: " <<
    setiosflags(ios::fixed) << setprecision(3) << "x[" <<
    gsr.axes.x(0) << ", " <<
    gsr.axes.x(1) << ", " <<
    gsr.axes.x(2) << "], y[" <<
    gsr.axes.y(0) << ", " <<
    gsr.axes.y(1) << ", " <<
    gsr.axes.y(2) << "], z[" <<
    gsr.axes.z(0) << ", " <<
    gsr.axes.z(1) << ", " <<
    gsr.axes.z(2) << "]");
  // initialize and save the Pt(A) curve
  LOG4CPLUS_DEBUG(logger, "initializing the Pt(A) curve");
  gsr.curve = GsrCurve(event, gsr.axes);
  gsr.APthCurve = Curve(gsr.curve.cols().x, event.dataNarrow().cols().Pth);
  LOG4CPLUS_DEBUG(logger, "initializing the branches of the Pt(A) curve");
  gsr.curve.initBranches("extremums").computeResidue();

  // middle time of the flux rope
  Time middleTime(event.config().beginTime);
  middleTime.add(event.config().endTime-event.config().beginTime, "second");
  // transform axes into Stonyhurst coordinate system
  if (event.config().spacecraft == "STA" ||
      event.config().spacecraft == "STB") {
    // convert from RTN to Stonyhurst
    Matrix3d rtn = getTransformationMatrix(event, "rtn", middleTime);
//    std::cout << rtn << std::endl<< std::endl;
    Matrix3d heeq = getTransformationMatrix(event, "heeq", middleTime);
    gsr.heeqAxes.x = heeq*(rtn.inverse()*gsr.axes.x);
    gsr.heeqAxes.y = heeq*(rtn.inverse()*gsr.axes.y);
    gsr.heeqAxes.z = heeq*(rtn.inverse()*gsr.axes.z);

  } else if (event.config().spacecraft == "WIND" ||
             event.config().spacecraft == "ACE") {
    // convert from GSE to Stonyhurst
    int retVal, es, es0;
    Vec zeros = {0,0,0};
    Vector3d heeq0 = Vector3d::Zero(3);
    es = date2es(middleTime.year(), middleTime.month(), middleTime.day(),
                 middleTime.hour(), middleTime.minute(), middleTime.second());
    es0 = date2es(middleTime.year(), middleTime.month(), middleTime.day(),
                  middleTime.hour(), middleTime.minute(), middleTime.second());
    retVal = cxform("GSE", "HEEQ", es0, zeros, heeq0.data());
    retVal = cxform("GSE", "HEEQ", es,
                    gsr.axes.x.data(), gsr.heeqAxes.x.data());
    gsr.heeqAxes.x -= heeq0;
    retVal = cxform("GSE", "HEEQ", es,
                    gsr.axes.y.data(), gsr.heeqAxes.y.data());
    gsr.heeqAxes.y -= heeq0;
    retVal = cxform("GSE", "HEEQ", es,
                    gsr.axes.z.data(), gsr.heeqAxes.z.data());
    gsr.heeqAxes.z -= heeq0;
  }

  LOG4CPLUS_DEBUG(logger, "HEEQ magnetic cloud axes: " <<
    setiosflags(ios::fixed) << setprecision(3) << "x[" <<
    gsr.heeqAxes.x(0) << ", " <<
    gsr.heeqAxes.x(1) << ", " <<
    gsr.heeqAxes.x(2) << "], y[" <<
    gsr.heeqAxes.y(0) << ", " <<
    gsr.heeqAxes.y(1) << ", " <<
    gsr.heeqAxes.y(2) << "], z[" <<
    gsr.heeqAxes.z(0) << ", " <<
    gsr.heeqAxes.z(1) << ", " <<
    gsr.heeqAxes.z(2) << "]");

  gsr.stonyhurstTheta = atan(gsr.heeqAxes.z(2)/
                        sqrt(pow(gsr.heeqAxes.z(0), 2)+
                             pow(gsr.heeqAxes.z(1), 2)))*180/M_PI;
  gsr.stonyhurstPhi = atan2(gsr.heeqAxes.z(1), gsr.heeqAxes.z(0))*180/M_PI;
  LOG4CPLUS_DEBUG(logger,
    "Invariant axis in Stonyhurst coordinates [theta, phi] = [" <<
    setiosflags(ios::fixed) << setprecision(1) <<
    gsr.stonyhurstTheta << ", " << gsr.stonyhurstPhi << "]");

  // calculate the magnetic field map and save it into gsr structure
  LOG4CPLUS_DEBUG(logger, "computing the magnetic field map");
  computeMap(event, gsr);

  // save the gsr results into event object
  event.gsr(gsr);
}

// run one loop of axes search algorithm
GsrResults GsrAnalyzer::loopAxes(Event& event,
                               double minTheta, double dTheta, double maxTheta,
                               double minPhi,   double dPhi,   double maxPhi)
{
  // get the logger instance
  Logger logger = Logger::getInstance("main");

  double theta, phi; // temporary theta and phi angles
  Axes axes; // temporary coordinate system
  GsrResults gsr; // holds the information about the current run

  // save axes limits in run structure
  gsr.minTheta = minTheta;
  gsr.dTheta   = dTheta;
  gsr.maxTheta = maxTheta;
  gsr.minPhi   = minPhi;
  gsr.dPhi     = dPhi;
  gsr.maxPhi   = maxPhi;

  // define residue matrix size
  gsr.originalResidue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                       int((maxPhi-minPhi)/dPhi)+1);
  gsr.combinedResidue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                       int((maxPhi-minPhi)/dPhi)+1);
  gsr.branchLength    = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                       int((maxPhi-minPhi)/dPhi)+1);

  MatrixXd pressureResidue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                            int((maxPhi-minPhi)/dPhi)+1);
  MatrixXd magneticResidue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                            int((maxPhi-minPhi)/dPhi)+1);

  // temporary quaternions for making axes rotations
  Quaterniond qTheta;
  Quaterniond qPhi;

  int i, k; // angle counters

  // temporary curve and axis required for the loop
  GsrCurve curve;
  Vector3d zTheta;

  i = 0; // starting from 1st row
  theta = minTheta; // initialize theta
  while (theta <= maxTheta) { // begin iteration through theta angles
    // initialize theta quaternion
    qTheta = AngleAxisd(theta*M_PI/180, event.pmvab().axes.y);
    zTheta = qTheta*event.pmvab().axes.z; // rotate z axis around PMVA y axis
    phi = minPhi; // initialize phi
    k = 0; // starting from 1st column
    while (phi <= maxPhi) { // begin iteration through phi angles
      // initialize phi quaternion
      qPhi = AngleAxisd(phi*M_PI/180, event.pmvab().axes.z);
      axes.z = (qPhi*zTheta).normalized(); // rotate z axis around PMVA z axis
      // initialize x axis
      axes.x = (event.dht().Vht.dot(axes.z)*axes.z-
                event.dht().Vht).normalized();
      // complement with y axis
      axes.y = axes.z.cross(axes.x).normalized();
      // create Pt(A) curve in new axes
      curve = GsrCurve(event, axes);

      Data data(event.dataNarrow());
      data.project(axes);
      BranchedCurve bzCurve(curve.cols().x, data.cols().Bz);
      bzCurve.initBranches("extremums").computeResidue();

      // initialize branches and compute the residue
      curve.initBranches("extremums").computeResidue();

      // save residue into a matrix
//      gsr.originalResidue(i,k) = bzCurve.originalResidue();//+curve.originalResidue();
//      gsr.combinedResidue(i,k) = bzCurve.combinedResidue();//+curve.combinedResidue();
//      gsr.branchLength(i,k) = bzCurve.branchLength();//+curve.branchLength();

      pressureResidue(i,k) = curve.combinedResidue();
      magneticResidue(i,k) = bzCurve.combinedResidue();

      gsr.originalResidue(i,k) = curve.originalResidue();
//      gsr.combinedResidue(i,k) = curve.combinedResidue();
      gsr.branchLength(i,k) = curve.branchLength();
      
      phi += dPhi; // make a step in phi
      k++; // move to the next column
    } // end iteration through phi angles
    theta += dTheta; // make a step in theta step
    i++; // move to the next row
  } // end iteration through theta angles

  gsr.combinedResidue = (pressureResidue.array()*magneticResidue.array()).matrix();

  return gsr; // return
}

// detect flux rope axes
GsrResults& GsrAnalyzer::detectAxes(Event& event, GsrResults& gsr) {

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  Data data(event.dataNarrow());
  if (data.symmetrize()) {
    // searching for y axis in symmetric case
    double sum, psum = 0, optYTheta = 0;
    for (double theta = -90; theta <= 90; theta+=0.1) {
      sum = 0;
      if (theta == 0) {
        sum = gsr.combinedResidue.col(0).array().inverse().sum()+
              gsr.combinedResidue.col(180).array().inverse().sum();
      } else {
        for (double phi = 0; phi <= 180; phi+=gsr.dPhi) {
          sum += 1/gsr.combinedResidue(
            floor(acos(cos(theta*M_PI/180)/
              sqrt(1+pow(sin(theta*M_PI/180)/tan((theta < 0 ? phi+180 : phi)*M_PI/180), 2)))*
              180/M_PI/gsr.dTheta+0.5),
            floor((theta < 0 ? phi+180 : phi)/gsr.dPhi+0.5));
        }
      }
      if (sum > psum) {
        psum = sum;
        optYTheta = theta;
      }
    }
    // quaternion, needed to turn to optimized Y axis
    Quaterniond qTheta;
    // initialize quaternion to get the Y axes of MC
    qTheta = AngleAxisd(-optYTheta*M_PI/180, event.pmvab().axes.x);
    // turn the axis
    gsr.axes.y = (qTheta*event.pmvab().axes.y).normalized();
    gsr.axes.x = -event.dht().Vht.normalized();
    gsr.axes.z = gsr.axes.x.cross(gsr.axes.y).normalized();

    cout << gsr.axes.y.dot(event.dht().Vht.normalized()) << endl;
    cout << gsr.axes.y.dot(Vector3d::UnitZ()) << endl;

    if (gsr.axes.y.dot(event.dht().Vht.normalized()) < 0.1) {
      cout << "HS" << endl;
      // quaternion, needed to cycle through all possible X axis orientations
      Quaterniond qPhi;
      // length of the SC intersection
      double L = (double)data.rows().size()/2;
      // SC intersection line length
      VectorXd x = VectorXd::LinSpaced(data.rows().size(), 0,
                                       data.rows().size()-1);
      x = (x.array()-L).matrix(); // center it, so it's symmetrical
      // project data on the trial axes
      data.project(gsr.axes);
      // temporary axes
      Axes axes = gsr.axes;

      // normalize Bx and Bz
      VectorXd Bx = (data.cols().Bx.array()/data.cols().By.array()*x.array()/L).matrix();
      VectorXd Bz = (data.cols().Bz.array()/data.cols().By.array()*x.array()/L).matrix();
      // remove spikes by applying median filtering
//      Bx = Filter::median1D(Bx, 25);
//      Bz = Filter::median1D(Bz, 25);
      Matrix2d M = Matrix2d::Zero(); // initialize 2x2 matrix with zeros
      M(0,0) = (Bz.array()*Bz.array()).matrix().mean()-Bz.mean()*Bz.mean();
      M(0,1) = (Bz.array()*Bx.array()).matrix().mean()-Bz.mean()*Bx.mean();
      M(1,0) = (Bx.array()*Bz.array()).matrix().mean()-Bx.mean()*Bz.mean();
      M(1,1) = (Bx.array()*Bx.array()).matrix().mean()-Bx.mean()*Bx.mean();
      // initialize eigen solver for adjoint matrix
      SelfAdjointEigenSolver<Matrix2d> eigensolver(M);
      // initialize and get eigen values in ascending order
      Vector2d eigenValues = eigensolver.eigenvalues();
      // initialize and get eigen vectors
      Matrix2d eigenVectors = eigensolver.eigenvectors();
      // maximum variance direction is the Z axis, minimum variance is X axis
      axes.x = (gsr.axes.z*eigenVectors.col(0)(0)+
                gsr.axes.x*eigenVectors.col(0)(1)).normalized();
      axes.z = (gsr.axes.z*eigenVectors.col(1)(0)+
                gsr.axes.x*eigenVectors.col(1)(1)).normalized();
      double criterionHS = eigenValues(1)/eigenValues(0);
      cout << eigenValues(0) << ", " << criterionHS << ", "
           << acos(axes.z.dot(event.pmvab().axes.z))*180/M_PI << ", "
           << atan2(axes.z.dot(event.pmvab().axes.y), axes.z.dot(event.pmvab().axes.x))*180/M_PI
           << endl;
      if (eigenValues(0) < 1e-3 and criterionHS > 3) {
        gsr.axes = axes;
      } else {
        cout << "MHS" << endl;
        // current minimum variance
        double lambda = 1e100;
        // current criterion
        double criterionMHS = 0;
        // iterate thtough kappa parameter [0,1]
        for (double k = 1; k >= 0; k-=0.001) {
          // normalize Bx and Bz
          VectorXd Bx = (data.cols().Bx.array()/data.cols().By.array()*(x.array()*k).sin()/k/L).matrix();
          VectorXd Bz = (data.cols().Bz.array()/data.cols().By.array()*(x.array()*k).sin()/k/L).matrix();
          // remove spikes by applying median filtering
//          Bx = Filter::median1D(Bx, 25);
//          Bz = Filter::median1D(Bz, 25);
          Matrix2d M = Matrix2d::Zero(); // initialize 2x2 matrix with zeros
          M(0,0) = (Bz.array()*Bz.array()).matrix().mean()-Bz.mean()*Bz.mean();
          M(0,1) = (Bz.array()*Bx.array()).matrix().mean()-Bz.mean()*Bx.mean();
          M(1,0) = (Bx.array()*Bz.array()).matrix().mean()-Bx.mean()*Bz.mean();
          M(1,1) = (Bx.array()*Bx.array()).matrix().mean()-Bx.mean()*Bx.mean();
          // initialize eigen solver for adjoint matrix
          SelfAdjointEigenSolver<Matrix2d> eigensolver(M);
          // initialize and get eigen values in ascending order
          Vector2d eigenValues = eigensolver.eigenvalues();
          if (eigenValues(0) <= lambda and eigenValues(1)/eigenValues(0) > criterionMHS) {
            lambda = eigenValues(0);
            criterionMHS = eigenValues(1)/eigenValues(0);
            // initialize and get eigen vectors
            Matrix2d eigenVectors = eigensolver.eigenvectors();
            // maximum variance direction is the Z axis, minimum variance is X axis
            axes.x = (gsr.axes.z*eigenVectors.col(0)(0)+
                      gsr.axes.x*eigenVectors.col(0)(1)).normalized();
            axes.z = (gsr.axes.z*eigenVectors.col(1)(0)+
                      gsr.axes.x*eigenVectors.col(1)(1)).normalized();
            cout << k << ", " << lambda << ", " << criterionMHS << ", "
                 << acos(axes.z.dot(event.pmvab().axes.z))*180/M_PI << ", "
                 << atan2(axes.z.dot(event.pmvab().axes.y), axes.z.dot(event.pmvab().axes.x))*180/M_PI
                 << endl;
          }
        }
        if (lambda < 1e-7/* && criterionMHS > criterionHS*/) {
          gsr.axes = axes;
        } else {
          goto residuals;
        }
      }
    }
  } else {
    residuals:
    cout << "GSR" << endl;
    // searching for the minimum residue direction
    int iTheta, iPhi; // index of optimal angles
    // get indices
    LOG4CPLUS_DEBUG(logger,
      "searching optimal axis in the range theta = [" <<
      event.config().minTheta << ", " << event.config().maxTheta << "], " <<
      "phi = [" << event.config().minPhi << ", " << event.config().maxPhi << "]");
    int iMinTheta = floor(event.config().minTheta/gsr.dTheta),
        iMaxTheta = ceil(event.config().maxTheta/gsr.dTheta)-1,
        iMinPhi   = floor(event.config().minPhi/gsr.dPhi),
        iMaxPhi   = ceil(event.config().maxPhi/gsr.dPhi)-1;

    gsr.combinedResidue.block(iMinTheta, iMinPhi,
                              iMaxTheta-iMinTheta+1, iMaxPhi-iMinPhi+1).
                        minCoeff(&iTheta, &iPhi);
    // translate indices into optimal angles
    gsr.optTheta = gsr.minTheta + (iMinTheta+iTheta)*gsr.dTheta;
    gsr.optPhi = gsr.minPhi + (iMinPhi+iPhi)*gsr.dPhi;

    // quaternions, needed to turn to optimized axes
    Quaterniond qTheta;
    Quaterniond qPhi;

    // initialize quaternions to turn the axes into MC axes
    qTheta = AngleAxisd(gsr.optTheta*M_PI/180, event.pmvab().axes.y);
    qPhi = AngleAxisd(gsr.optPhi*M_PI/180, event.pmvab().axes.z);
    // turn the axes
    gsr.axes.z = (qPhi*(qTheta*event.pmvab().axes.z)).normalized();
    gsr.axes.x = (event.dht().Vht.dot(gsr.axes.z)*gsr.axes.z-
                  event.dht().Vht).normalized();
    gsr.axes.y = gsr.axes.z.cross(gsr.axes.x).normalized();
    
//    std::cout << gsr.axes.x << std::endl << std::endl;
//    std::cout << gsr.axes.y << std::endl << std::endl;
//    std::cout << gsr.axes.z << std::endl << std::endl;
  }
//  std::cout << event.dht().Vht.dot(gsr.axes.x) << std::endl << std::endl;
/*
  if (event.dht().Vht.dot(gsr.axes.x) > 0) {
    gsr.axes.z = -gsr.axes.z;
    gsr.axes.x = -gsr.axes.x;
  }

  if (acos(gsr.axes.z.dot(event.pmvab().axes.z)) > M_PI/2) {
    gsr.axes.y = -gsr.axes.y;
    gsr.axes.z = -gsr.axes.z;
  }
*/
  gsr.optTheta = acos(gsr.axes.z.dot(event.pmvab().axes.z))*180/M_PI;
  gsr.optPhi = atan2(gsr.axes.z.dot(event.pmvab().axes.y), gsr.axes.z.dot(event.pmvab().axes.x))*180/M_PI;

  LOG4CPLUS_DEBUG(logger,
      "optimal angles for the invariant axes [theta, phi] = [" <<
      gsr.optTheta << ", " << gsr.optPhi << "]");

  return gsr;

}

// detect flux rope axes in the non-symmetrical case
GsrResults& GsrAnalyzer::detectAxesByResidue(Event& event, GsrResults& gsr) {

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  // searching for the minimum residue direction
  int iTheta, iPhi; // index of optimal angles
  // get indices
  LOG4CPLUS_DEBUG(logger,
    "searching optimal axis in the range theta = [" <<
    event.config().minTheta << ", " << event.config().maxTheta << "], " <<
    "phi = [" << event.config().minPhi << ", " << event.config().maxPhi << "]");
  int iMinTheta = floor(event.config().minTheta/gsr.dTheta),
      iMaxTheta = ceil(event.config().maxTheta/gsr.dTheta)-1,
      iMinPhi   = floor(event.config().minPhi/gsr.dPhi),
      iMaxPhi   = ceil(event.config().maxPhi/gsr.dPhi)-1;

  gsr.combinedResidue.block(iMinTheta, iMinPhi,
                            iMaxTheta-iMinTheta+1, iMaxPhi-iMinPhi+1).
                      minCoeff(&iTheta, &iPhi);
  // translate indices into optimal angles
  gsr.optTheta = gsr.minTheta + (iMinTheta+iTheta)*gsr.dTheta;
  gsr.optPhi = gsr.minPhi + (iMinPhi+iPhi)*gsr.dPhi;

  LOG4CPLUS_DEBUG(logger,
    "optimal angles for the invariant axes [theta, phi] = [" <<
    gsr.optTheta << ", " << gsr.optPhi << "]");

  // quaternions, needed to turn to optimized axes
  Quaterniond qTheta;
  Quaterniond qPhi;

  // initialize quaternions to turn the axes into MC axes
  qTheta = AngleAxisd(gsr.optTheta*M_PI/180, event.pmvab().axes.y);
  qPhi = AngleAxisd(gsr.optPhi*M_PI/180, event.pmvab().axes.z);
  // turn the axes
  gsr.axes.z = (qPhi*(qTheta*event.pmvab().axes.z)).normalized();
  gsr.axes.x = (event.dht().Vht.dot(gsr.axes.z)*gsr.axes.z-
                event.dht().Vht).normalized();
  gsr.axes.y = gsr.axes.z.cross(gsr.axes.x).normalized();

  return gsr; // return
}

// detect flux rope axes in case of symmetry
GsrResults& GsrAnalyzer::detectAxesByVariance(Event& event, GsrResults& gsr) {

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  // searching for y axis in symmetric case
  double sum, psum = 0, optYTheta = 0;
  for (double theta = -90; theta <= 90; theta+=0.1) {
    sum = 0;
    if (theta == 0) {
      sum = gsr.combinedResidue.col(0).array().inverse().sum()+
            gsr.combinedResidue.col(180).array().inverse().sum();
    } else {
      for (double phi = 0; phi <= 180; phi+=gsr.dPhi) {
        sum += 1/gsr.combinedResidue(
          floor(acos(cos(theta*M_PI/180)/
            sqrt(1+pow(sin(theta*M_PI/180)/tan((theta < 0 ? phi+180 : phi)*M_PI/180), 2)))*
            180/M_PI/gsr.dTheta+0.5),
          floor((theta < 0 ? phi+180 : phi)/gsr.dPhi+0.5));
      }
    }
    if (sum > psum) {
      psum = sum;
      optYTheta = theta;
    }
  }

  // quaternion, needed to turn to optimized Y axis
  Quaterniond qTheta;
  // initialize quaternion to get the Y axes of MC
  qTheta = AngleAxisd(-optYTheta*M_PI/180, event.pmvab().axes.x);
  // turn the axis
  gsr.axes.y = (qTheta*event.pmvab().axes.y).normalized();
  gsr.axes.x = -event.dht().Vht.normalized();
  gsr.axes.z = gsr.axes.x.cross(gsr.axes.y).normalized();

  // we iterate through all possible orientations of X axis and estimate
  // variances
  // quaternion, needed to cycle through all possible X axis orientations
  Quaterniond qPhi;
  // initialize the Pt(A) curve and it's branches to get the indices of the
  // symmetric interval
  GsrCurve curve = GsrCurve(event, gsr.axes); // create Pt(A) curve
  curve.initBranches("extremums"); // find it's branches
  // length of the SC intersection
  double L = (double)(curve.rightIndex()-curve.leftIndex())/2;
  // SC intersection line length
  VectorXd x = VectorXd::LinSpaced(curve.rightIndex()-curve.leftIndex()+1, 0,
                                   curve.rightIndex()-curve.leftIndex());
  x = (x.array()-L).matrix(); // center it, so it's symmetrical

  Data data(event.dataNarrow());
  // project data on the trial axes
  data.project(gsr.axes);
  // current minimum variance
  double lambda = 1e100;
  // filter data to get symmetric interval
  data.filter(curve.leftIndex(), curve.rightIndex());
  // temporary axes
  Axes axes = gsr.axes;

//  Bx.segment(curve.centerIndex()-10-1,
//             curve.rightIndex()-curve.centerIndex()) =
//  Bx.segment(curve.centerIndex()+10+1,
//             curve.rightIndex()-curve.centerIndex());

//  Bz.segment(curve.centerIndex()-10-1,
//             curve.rightIndex()-curve.centerIndex()) =
//  Bz.segment(curve.centerIndex()+10+1,
//             curve.rightIndex()-curve.centerIndex());

  // remove spikes by applying median filtering
//  Bx = Filter::median1D(Bx, 25);
//  Bz = Filter::median1D(Bz, 25);

  // iterate thtough kappa parameter [0,1]
  for (double k = 1; k >= 0; k-=0.001) {
    // normalize Bx and Bz
    VectorXd Bx = (data.cols().Bx.array()/data.cols().By.array()*(x.array()*k).sin()/k/L).matrix();
    VectorXd Bz = (data.cols().Bz.array()/data.cols().By.array()*(x.array()*k).sin()/k/L).matrix();

//    int n = 1;

//    Bx.segment(curve.centerIndex()-n,
//               curve.rightIndex()-curve.centerIndex()-n) =
//    Bx.tail(curve.rightIndex()-curve.centerIndex()-n);
//    Bx.conservativeResize(Bx.size()-1-n*2);

//    Bz.segment(curve.centerIndex()-n,
//               curve.rightIndex()-curve.centerIndex()-n) =
//    Bz.tail(curve.rightIndex()-curve.centerIndex()-n);
//    Bz.conservativeResize(Bz.size()-1-n*2);

    // remove spikes by applying median filtering
    Bx = Filter::median1D(Bx, 25);
    Bz = Filter::median1D(Bz, 25);

    Matrix2d M = Matrix2d::Zero(); // initialize 2x2 matrix with zeros
    M(0,0) = (Bz.array()*Bz.array()).matrix().mean()-Bz.mean()*Bz.mean();
    M(0,1) = (Bz.array()*Bx.array()).matrix().mean()-Bz.mean()*Bx.mean();
    M(1,0) = (Bx.array()*Bz.array()).matrix().mean()-Bx.mean()*Bz.mean();
    M(1,1) = (Bx.array()*Bx.array()).matrix().mean()-Bx.mean()*Bx.mean();
//    cout << M << endl;
    // initialize eigen solver for adjoint matrix
    SelfAdjointEigenSolver<Matrix2d> eigensolver(M);
    // initialize and get eigen values in ascending order
    Vector2d eigenValues = eigensolver.eigenvalues();
//    cout << eigenValues(0) << " -> " << eigenValues(1) << endl;
    if (eigenValues(0) <= lambda) {
      lambda = eigenValues(0);
      // initialize and get eigen vectors
      Matrix2d eigenVectors = eigensolver.eigenvectors();
      // maximum variance direction is the Z axis, minimum variance is X axis
      axes.x = (gsr.axes.z*eigenVectors.col(0)(0)+
                gsr.axes.x*eigenVectors.col(0)(1)).normalized();
      axes.z = (gsr.axes.z*eigenVectors.col(1)(0)+
                gsr.axes.x*eigenVectors.col(1)(1)).normalized();
      cout << k << ", " << lambda << ", "
           << acos(axes.z.dot(event.pmvab().axes.z))*180/M_PI << ", "
           << endl;
//      Plotter plotter;
//      plotter.plotData1D(Bx);
//      plotter.plotData1D(Bz);
    }
  }
  gsr.axes = axes;
  /*
  double var = 0, pvar = 1e9; // current and previous variance
  // iterate through rotation angles
  for (double phi = -90; phi <= 90; phi+=0.1) {
    // initalize the quaternion
    qPhi = AngleAxisd(phi*M_PI/180, gsr.axes.y);
    // rotate X axis to get a new trial one
    axes.x = (qPhi*gsr.axes.x).normalized();
    // calculate new trial Z axis
    axes.z = axes.x.cross(axes.y).normalized();
    // geta copy of data
    Data data(event.dataNarrow());
    // project data on the trial axes
    data.project(axes);
    // filter data to get symmetric interval
    data.filter(curve.leftIndex(), curve.rightIndex());
    // normalize Bx
    VectorXd Bx = (data.cols().Bx.array()/data.cols().By.array()*
                   x.array()/L).matrix();
    // remove spikes by applying median filtering
    Bx = Filter::median1D(Bx, 25);
    // estimate variance of Bx
    double var = (Bx.array()-Bx.mean()).pow(2).mean();

    if (var < pvar) {
      pvar = var;
      faxes = axes;
    }
  }
  gsr.axes = faxes;
  */

  if (event.dht().Vht.dot(gsr.axes.x) > 0) {
    gsr.axes.z = -gsr.axes.z;
    gsr.axes.x = -gsr.axes.x;
  }

  if (acos(gsr.axes.z.dot(event.pmvab().axes.z)) > M_PI/2) {
    gsr.axes.y = -gsr.axes.y;
    gsr.axes.z = -gsr.axes.z;
  }

  gsr.optTheta = acos(gsr.axes.z.dot(event.pmvab().axes.z))*180/M_PI;
  gsr.optPhi = atan2(gsr.axes.z.dot(event.pmvab().axes.y), gsr.axes.z.dot(event.pmvab().axes.x))*180/M_PI;

  LOG4CPLUS_DEBUG(logger,
    "optimal angles for the invariant axes [theta, phi] = [" <<
    gsr.optTheta << ", " << gsr.optPhi << "]");

  return gsr;
}

// compute magnetic field map for a given GSR run
GsrResults& GsrAnalyzer::computeMap(Event& event, GsrResults& gsr) {

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  // copy the data into external Dara object, we need a copy because we
  // will have to project it onto new coordinate system
  Data data(event.dataNarrow());

//  std::cout << Curve::resampled(data.cols().Bx, 15) << std::endl << std::endl;
//  std::cout << Curve::resampled(data.cols().By, 15) << std::endl << std::endl;
//  std::cout << Curve::resampled(data.cols().Bz, 15) << std::endl << std::endl;

  // project the data into optimized coordinate system
  data.project(gsr.axes);

//  std::cout << Curve::resampled(data.cols().Bx, 15) << std::endl << std::endl;
//  std::cout << Curve::resampled(data.cols().By, 15) << std::endl << std::endl;
//  std::cout << Curve::resampled(data.cols().Bz, 15) << std::endl << std::endl;

//  cout << data.cols().Vx/1e3 << endl << endl;
//  cout << data.cols().Vy/1e3 << endl << endl;

  // compute the dx step
  double dx = -event.dht().Vht.dot(gsr.axes.x)*event.config().samplingInterval*
               gsr.curve.size()/event.config().Nx;
  // first estimate for the dy step
  double dy = event.config().ratio*dx;
  // number of Y-nodes
  int Ny = floor((event.config().maxY-event.config().minY)*
           GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy);
  // final value for the Y step
  dy = (event.config().maxY-event.config().minY)*
       GSL_CONST_MKSA_ASTRONOMICAL_UNIT/Ny;
  // number of steps into positive Y direction
  int NyUp = ceil(event.config().maxY*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy);
  // number of steps into negative Y direction
  int NyDown = -floor(event.config().minY*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy);
  // final number of Y-nodes
  Ny = NyDown + 1 + NyUp;
  // initialize the vector of Y nodes
  VectorXd Y = VectorXd::Zero(Ny);
  // and fill it with actual Y values
  for (int i = -NyDown; i <= NyUp; i++) {
    Y(NyDown+i) = i*dy;
  }
  // number of X nodes
  int Nx = event.config().Nx;

  // initialize and fill in the X coordinates vector (on a lattice)
  VectorXd X = VectorXd::Zero(Nx);
  for (int i = 0; i < Nx; i++) {
    X(i) = i*dx;
  }

  // save lattice parameters in the run structure
  gsr.dx = dx;
  gsr.dy = dy;
  gsr.Nx = Nx;
  gsr.Ny = Ny;
  gsr.X = X;
  gsr.Y = Y;

  // copy branches of Pt(A) into separate Curve objects, we need the copies
  // because we will need to resample them
  Curve APtInCurve(gsr.curve.branches()[0]);
  Curve APtOutCurve(gsr.curve.branches()[1]);

  // resample the branches curves to the new lattice
  APtInCurve.resample(Nx);
  APtOutCurve.resample(Nx);

  // define Bx and By for quiver plot
  gsr.Bx = Curve::resampled(data.cols().Bx, Nx);
//  std::cout << gsr.Bx << std::endl << std::endl;
  gsr.By = Curve::resampled(data.cols().By, Nx);
//  std::cout << gsr.By << std::endl << std::endl;

  // initialize arrays for the branches
  double *AarrAll  = new double[2*Nx],
          *PtArrAll = new double[2*Nx];
  // fill the arrays with values of the both branches
  for (int i = 0; i < 2*Nx; i++) {
    if (i < Nx) { // get values from the inward branch
      AarrAll[i] = APtInCurve.cols().x(i);
      PtArrAll[i] = APtInCurve.cols().y(i);
    } else { // get values from the outward branch
      AarrAll[i] = APtOutCurve.cols().x(i-Nx);
      PtArrAll[i] = APtOutCurve.cols().y(i-Nx);
    }
  }

  // map arrays into vectors
  Map<VectorXd> Aall(AarrAll, 2*Nx);
  Map<VectorXd> PtAll(PtArrAll, 2*Nx);

  // build a Pt(A) curve from the data vectors
  Curve APtAllCurve(Aall, PtAll);

  // sort and unique the curve
  APtAllCurve.sort().unique();

  int nAll = APtAllCurve.size(); // length of the curve

  // change size of the data arrays
  AarrAll = (double*) realloc(AarrAll, nAll*sizeof(double));
  PtArrAll = (double*) realloc(PtArrAll, nAll*sizeof(double));
  // reinitialize the data arrays
  AarrAll  = const_cast<double*>(APtAllCurve.cols().x.data());
  PtArrAll = const_cast<double*>(APtAllCurve.cols().y.data());

  // polyexp fitting
  PolyExpFit APtFit(nAll, AarrAll, PtArrAll, event.config().order,
    event.config().fittingParameterCtr, event.config().fittingParameterBdr);
  APtFit.fit();

  // find the slope of the curve
  PolyFit APtLineFit(nAll, AarrAll, PtArrAll, 1);
  APtLineFit.fit();
  // initialize the slope, TODO: sign function
  int slope = (APtLineFit.c()[1] > 0 ? 1 : -1);

  // do the fitting
  PolyFit APtPolyFit(nAll, AarrAll, PtArrAll, event.config().order);
  APtPolyFit.fit();

  double Ab, Ac; // boundary and center value of A
  VectorXd Atmp; // wide range A, just for plotting
  int AbIndex; // index of the
  if (slope > 0) { // Ab < Ac
    Ab = APtAllCurve.cols().x.minCoeff(&AbIndex);
    Ac = APtAllCurve.cols().x.maxCoeff();
  } else { // Ac < Ab
    Ab = APtAllCurve.cols().x.maxCoeff(&AbIndex);
    Ac = APtAllCurve.cols().x.minCoeff();
  }

  // do not allow the Ab to pass the zero derivative point of the polynomial
  while (slope > 0 && APtPolyFit.df(Ab) < 0 ||
         slope < 0 && APtPolyFit.df(Ab) > 0)
  {
    if (slope > 0) {
      Ab = APtAllCurve.cols().x(++AbIndex);
    } else {
      Ab = APtAllCurve.cols().x(--AbIndex);
    }
  }

//  Ab = 53;

  // determine limits for Pt(A) plots
  double AcLim,
         AcLim1,
         AcLim2,
         PtMax = 1.5*APtFit.f(Ac),
         deltaPt = 0.01*(APtFit.f(Ac)-APtFit.f(Ab));
  if (slope > 0) {
    AcLim1 = Ac;
    AcLim2 = Ac+(Ac-Ab);
    AcLim  = AcLim2;
  } else {
    AcLim1 = Ac-(Ab-Ac);
    AcLim2 = Ac;
    AcLim  = AcLim1;
  }
  if (APtFit.f(AcLim) > PtMax) {
    while (abs(APtFit.f(AcLim)-PtMax) > deltaPt) {
      AcLim = (AcLim1+AcLim2)/2;
      if (APtFit.f(AcLim) > PtMax) {
        if (slope > 0) {
          AcLim2 = AcLim;
        } else {
          AcLim1 = AcLim;
        }
      } else {
        if (slope > 0) {
          AcLim1 = AcLim;
        } else {
          AcLim2 = AcLim;
        }
      }
    }
  }
  double AbLim = (slope > 0 ? Ab-3*(Ac-Ab) : Ab+3*(Ab-Ac));
  if (slope > 0) {
    Atmp = VectorXd::LinSpaced(1000, AbLim, AcLim);
  } else {
    Atmp = VectorXd::LinSpaced(1000, AcLim, AbLim);
  }

  // save central and boundary values of the vector potential
  gsr.Ac = Ac;
  gsr.Ab = Ab;
///*!!!*/  gsr.Ab = Ac;
  gsr.beginTime = event.config().beginTime;
  gsr.endTime = event.config().endTime;

  for (int i = 0; i < gsr.curve.cols().x.size(); i++) {
    if (gsr.Ac > gsr.Ab && gsr.curve.cols().x(i) > gsr.Ab ||
        gsr.Ac < gsr.Ab && gsr.curve.cols().x(i) < gsr.Ab) {
      gsr.beginTime = Time(event.dataNarrow().row(i).year,
                           event.dataNarrow().row(i).month,
                           event.dataNarrow().row(i).day,
                           event.dataNarrow().row(i).hour,
                           event.dataNarrow().row(i).minute,
                           event.dataNarrow().row(i).second);
      break;
    }
  }
  for (int i = gsr.curve.cols().x.size()-1; i >= 0; i--) {
    if (gsr.Ac > gsr.Ab && gsr.curve.cols().x(i) > gsr.Ab ||
        gsr.Ac < gsr.Ab && gsr.curve.cols().x(i) < gsr.Ab) {
      gsr.endTime = Time(event.dataNarrow().row(i).year,
                         event.dataNarrow().row(i).month,
                         event.dataNarrow().row(i).day,
                         event.dataNarrow().row(i).hour,
                         event.dataNarrow().row(i).minute,
                         event.dataNarrow().row(i).second);
      break;
    }
  }
  LOG4CPLUS_DEBUG(logger, "GSR boundaries: " << setfill('0') <<
    gsr.beginTime.year() << '-' <<
    setw(2) << gsr.beginTime.month() << '-' <<
    setw(2) << gsr.beginTime.day() << ' ' <<
    setw(2) << gsr.beginTime.hour() << ':' <<
    setw(2) << gsr.beginTime.minute() << ':' <<
    setw(2) << gsr.beginTime.second() << " - " <<
    gsr.endTime.year() << '-' <<
    setw(2) << gsr.endTime.month() << '-' <<
    setw(2) << gsr.endTime.day() << ' ' <<
    setw(2) << gsr.endTime.hour() << ':' <<
    setw(2) << gsr.endTime.minute() << ':' <<
    setw(2) << gsr.endTime.second());

  // initialize temporary Pt vector, for plotting only
  VectorXd PtTmp = VectorXd::Zero(1000);
  // fill it with evaluated values for the fitted Pt(A)
  for (int i = 0; i < 1000; i++) {
    // polyexp fitting
    PtTmp(i) = APtFit.f(Atmp[i]);
  }

  // initialize temporary fitted curve, for plotting only
  Curve APtFitCurve(Atmp, PtTmp);

  gsr.APtInCurve  = APtInCurve;
  gsr.APtOutCurve = APtOutCurve;
  gsr.APtFitCurve = APtFitCurve;

  // initialize temporary dPt/dA vector, for plotting only
  VectorXd dPtTmp = VectorXd::Zero(1000);
  // fill it with derivatives
  for (int i = 0; i < 1000; i++) {
    // polyexp fitting
    dPtTmp(i) = APtFit.df(Atmp[i]);
  }

  // initialize temporary derivative of the fitted curve, for plotting only
  Curve AdPtFitCurve(Atmp, dPtTmp);

  gsr.AdPtFitCurve = AdPtFitCurve;

  // copy A, B and Pth to new vector objects, needed for resampling
  VectorXd A(gsr.curve.cols().x);
  VectorXd Bx(data.cols().Bx);
  VectorXd By(data.cols().By);
  VectorXd Bz(data.cols().Bz);
  VectorXd Pth(data.cols().Pth);

  // resample physical parameters to use on a lattice
  Curve::resample(A, Nx);
  Curve::resample(Bx, Nx);
  Curve::resample(By, Nx);
  Curve::resample(Bz, Nx);
  Curve::resample(Pth, Nx);

  // initialize matrices for A, Bx and Bz maps
  MatrixXd Axy = MatrixXd::Zero(Ny, Nx),
           Bxy = MatrixXd::Zero(Ny, Nx),
           Bzy = MatrixXd::Zero(Ny, Nx);

  // fill initialize A and Bx values to the matrices, i.e. measured spacecraft
  // data
  Axy.row(NyDown) = A;
//  std::cout << Bx << std::endl << std::endl;
  Bxy.row(NyDown) = Bx;

  // initialize vectors for some important derivatives
  VectorXd d2A_dx2 = VectorXd::Zero(Nx),
           dPt_dA  = VectorXd::Zero(Nx),
           d2A_dy2 = VectorXd::Zero(Nx);

  // initialize differentiator
  Differentiator differentiator;

  LOG4CPLUS_DEBUG(logger, NyUp << " steps up and " << NyDown <<
                          " steps down with step size " << dy);
  // reconstruct the upper part of the map using recursive solver
  LOG4CPLUS_DEBUG(logger, "reconstructing the upper part of the map");
  for (int i = 1; i <= NyUp; i++) {
    // 2nd derivative of A by x using Holoborodko2 filter, smoothed with
    // weighted average prior to differenting
    d2A_dx2 = differentiator.Holoborodko2(5,Axy.row(NyDown+i-1),dx);
    // evaluate the 1st derivative of Pt by A using fitting curve
    for (int k = 0; k < Nx; k++) {
      // polyexp fitting
      dPt_dA(k) = APtFit.df(Axy(NyDown+i-1, k));
    }
    // compute 2nd derivative of A by y using Grad-Shafranov equation
    d2A_dy2 = -d2A_dx2-GSL_CONST_MKSA_VACUUM_PERMEABILITY*dPt_dA;
    // write the next row of A, GS equation used
    Axy.row(NyDown+i) = Axy.row(NyDown+i-1)+Bxy.row(NyDown+i-1)*dy+
                        d2A_dy2.transpose()*pow(dy, 2)/2;
    // smooth it with weighted average
    Axy.row(NyDown+i) = Curve::weightedAverage(Axy.row(NyDown+i),
                                               min(double(i)/double(NyUp),0.7));
    // wtite the next row of Bx
    Bxy.row(NyDown+i) = Bxy.row(NyDown+i-1)+d2A_dy2.transpose()*dy;
    // smooth it with weighted average
    Bxy.row(NyDown+i) = Curve::weightedAverage(Bxy.row(NyDown+i),
                                               min(double(i)/double(NyUp),0.7));
//    std::cout << Bxy.row(NyDown+i) << std::endl << std::endl;
  }

  // reconstruct the lower part of teh map using recursive solver,
  // everything is the same as for the upper part
  LOG4CPLUS_DEBUG(logger, "reconstructing the lower part of the map");
  for (int i = -1; i >= -NyDown; i--) {
    d2A_dx2 = differentiator.Holoborodko2(5,Axy.row(NyDown+i+1),dx);
    for (int k = 0; k < Nx; k++) {
      // polyexp fitting
      dPt_dA(k) = APtFit.df(Axy(NyDown+i+1, k));
    }
    d2A_dy2 = -d2A_dx2-GSL_CONST_MKSA_VACUUM_PERMEABILITY*dPt_dA;
    Axy.row(NyDown+i) = Axy.row(NyDown+i+1)-Bxy.row(NyDown+i+1)*dy+
                        d2A_dy2.transpose()*pow(dy, 2)/2;
    Axy.row(NyDown+i) = Curve::weightedAverage(Axy.row(NyDown+i),
                                               min(double(i)/double(-NyDown),0.7));
    Bxy.row(NyDown+i) = Bxy.row(NyDown+i+1)-d2A_dy2.transpose()*dy;
    Bxy.row(NyDown+i) = Curve::weightedAverage(Bxy.row(NyDown+i),
                                               min(double(i)/double(-NyDown),0.7));
  }

  gsr.Axy = Axy; // save vector potential

  int iRow, iCol;
  if (slope > 0) {
    gsr.Aa = gsr.Axy.maxCoeff(&iRow, &iCol);
  } else {
    gsr.Aa = gsr.Axy.minCoeff(&iRow, &iCol);
  }
  gsr.ip = gsr.Y(iRow)*copysign(1.0, gsr.axes.y(2));

  LOG4CPLUS_DEBUG(logger, "Impact parameter = " <<
    gsr.ip/GSL_CONST_MKSA_ASTRONOMICAL_UNIT << "AU");

  LOG4CPLUS_DEBUG(logger, "estimating Bz all over the map");

//  (gsr.Axy.row(0).head(iCol+1)-gsr.Ab).abs().minCoeff(&iCol)
//  (gsr.Axy.row(0).head(iCol+1)-gsr.Ab).abs().minCoeff(&iCol)

  // initialize the Bz(A) curve
  Curve ABzCurve(A, Bz);

  // fit it with exponent
  ExpFit ABzFit(Nx, A.data(), Bz.data());
//  PosZeroExpFit ABzFit(Nx, A.data(), Bz.data());
//  PolyFit ABzFit(Nx, A.data(), Bz.data(), 2);
//  PolyExpFit ABzFit(Nx, A.data(), Bz.data(), 2, 2, 2);

//  Curve ABzCurveCopy = ABzCurve;
//  ABzCurveCopy.sort().unique();
//  PolyExpFit ABzFit(Nx, const_cast<double*>(ABzCurveCopy.cols().x.data()),
//                    const_cast<double*>(ABzCurveCopy.cols().y.data()),
//                    event.config().order,
//                    event.config().fittingParameterCtr,
//                    event.config().fittingParameterBdr);

  ABzFit.fit();
  VectorXd BzFit = VectorXd::Zero(Nx); // vector of fitted Bz values
  // fill it by evaluating the fitting function
  for (int i = 0; i < Nx; i++) {
    BzFit(i) = ABzFit.f(A(i));
  }

  // initialize fitted Bz(A)
  Curve ABzFitCurve(A, BzFit);
  // sort and unique it
  ABzFitCurve.sort().unique();

  // save Bz(A)
  gsr.ABzCurve = ABzCurve;
  gsr.ABzFitCurve = ABzFitCurve;

  // calculate the Bz map using the fitted Bz(A)
  for (int i = -NyDown; i <= NyUp; i++) {
    for (int k = 0; k < Nx; k++) {
      Bzy(NyDown+i, k) = ABzFit.f(Axy(NyDown+i, k));
    }
  }

  gsr.Bz = Bzy; // save magnetic field map

  // free dynamic arrays
//  delete [] AarrAll;
//  delete [] PtArrAll;
}

// read transformation matrix
Matrix3d GsrAnalyzer::getTransformationMatrix(Event& event, string cs,
                                              Time middleTime)
{
  ifstream dataFileStream; // stream from the file with data
  string dataFileLine; // a single line from a file as a string
  istringstream dataFileLineStream; // a stream from a line from a file
  ostringstream dataPathStream; // a stream for a data file path
  Time currentTime; // Time object for storing current time for comparison
  int delta = 5*24*3600;
  int year, doy, second, flag;

  Matrix3d tm;

  // determine the path to the data file dependant on spacecraft
  // push the path to data files
  dataPathStream << event.dataDir() << '/';
  if (event.config().spacecraft == "STA") {
    dataPathStream << "stereo_a_";
  } else if (event.config().spacecraft == "STB") {
    dataPathStream << "stereo_b_";
  } else { // throw an error - unknown spacecraft

  }
  boost::algorithm::to_lower(cs);
  dataPathStream << cs << ".dat";
  // stringstream to string
  string dataPath = dataPathStream.str();
  // open data file as a stream
  dataFileStream.open(dataPath.c_str());
  // check if the file was opened correctly
  if (dataFileStream.is_open()) {
//    std::cout << "CS file opened" << std::endl;
    // start iterating through data file line, one timestamp at a time
    while (dataFileStream.good()) {
      // get line from the file as a string
      getline(dataFileStream, dataFileLine);
      // check if the line contains actually data
      if (!dataFileLine.empty() && // check that it's not empty
          dataFileLine[0] != '#') // and not commented out
      {
        // clear string stream from the previus line of data
        dataFileLineStream.clear();
        // use the next line of data as a source for the string stream
        dataFileLineStream.str(dataFileLine);
        // parse the data from the stream of line of the data file
        dataFileLineStream >> year >> doy >> second >> flag;
        // initialize current Time object with time data
//        std::cout << year << "--" << second << std::endl;
        currentTime = Time(year, doy, 0, 0, 0.0);
        currentTime.add(second, "second");
//        std:: cout << abs(currentTime-middleTime) << std::endl;
        if (abs(currentTime-middleTime) < delta) {
//          std::cout << "CS delta" << std::endl;
          delta = abs(currentTime-middleTime);
          dataFileLineStream >>
            tm(0,0) >> tm(0,1) >> tm(0,2) >>
            tm(1,0) >> tm(1,1) >> tm(1,2) >>
            tm(2,0) >> tm(2,1) >> tm(2,2);
//            std::cout << "-- " <<  tm(0,0) << std::endl;
        } else {
          continue;
        }
      }
    } // end of iteration through the lines of the data file
  }
  return tm;
}


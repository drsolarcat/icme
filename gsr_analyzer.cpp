
// project headers
#include "gsr_analyzer.h"
#include "dht_analyzer.h"
#include "mva_analyzer.h"
#include "gsr_curve.h"
#include "fit_poly.h"
#include "fit_exp.h"
#include "fit_poly_exp.h"
#include "differentiator.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_fit.h>
#include <log4cplus/logger.h>
#include <log4cplus/configurator.h>
// standard headers
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;
using namespace log4cplus;

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

  // make a run of axes searching algorithm, save the results in the gsr
  // structure
  LOG4CPLUS_DEBUG(logger, "searching for optimal axes with 1 degree step");
  GsrResults gsr = loopAxes(event, 0, 1, 90, 0, 1, 360);

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
  // flip axes if necessary
  Data dataTmp(event.dataNarrow());
  dataTmp.project(gsr.axes);
  VectorXd xTmp = VectorXd::LinSpaced(dataTmp.rows().size(), 0,
                                      dataTmp.rows().size()-1);
  PolyFit BzTmpFit(dataTmp.rows().size(), xTmp.data(),
                   const_cast<double*>(dataTmp.cols().Bz.data()), 2);
  BzTmpFit.fit();
  if (BzTmpFit.c()[2] > 0) {
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
  LOG4CPLUS_DEBUG(logger, "initializing the branches of the Pt(A) curve");
  gsr.curve.initBranches("extremums").computeResidue();

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
      // initialize branches and compute the residue
      curve.initBranches("extremums").computeResidue();
      // save residue into a matrix
      gsr.originalResidue(i,k) = curve.originalResidue();
      gsr.combinedResidue(i,k) = curve.combinedResidue();
      gsr.branchLength(i,k) = curve.branchLength();
      phi += dPhi; // make a step in phi
      k++; // move to the next column
    } // end iteration through phi angles
    theta += dTheta; // make a step in theta step
    i++; // move to the next row
  } // end iteration through theta angles

  // searching for the minimum residue direction
  int iTheta, iPhi; // index of optimal angles
  // get indices
  LOG4CPLUS_DEBUG(logger,
    "searching optimal axis in the range theta = [" <<
    event.config().minTheta << ", " << event.config().maxTheta << "], " << 
    "phi = [" << event.config().minPhi << ", " << event.config().maxPhi << "]");
  int iMinTheta = floor(event.config().minTheta/dTheta),
      iMaxTheta = ceil(event.config().maxTheta/dTheta)-1,
      iMinPhi   = floor(event.config().minPhi/dPhi),
      iMaxPhi   = ceil(event.config().maxPhi/dPhi)-1;
  gsr.combinedResidue.block(iMinTheta, iMinPhi, 
                            iMaxTheta-iMinTheta+1, iMaxPhi-iMinPhi+1).
                      minCoeff(&iTheta, &iPhi);
  // translate indices into optimal angles
  gsr.optTheta = minTheta + (iMinTheta+iTheta)*dTheta;
  gsr.optPhi = minPhi + (iMinPhi+iPhi)*dPhi;

  LOG4CPLUS_DEBUG(logger,
    "optimal angles for the invariant axes [theta, phi] = [" <<
    gsr.optTheta << ", " << gsr.optPhi << "]");

  return gsr; // return
}

// compute magnetic field map for a given GSR run
GsrResults& GsrAnalyzer::computeMap(Event& event, GsrResults& gsr) {

  // get the logger instance
  Logger logger = Logger::getInstance("main");

  // copy the data into external Dara object, we need a copy because we
  // will have to project it onto new coordinate system
  Data data(event.dataNarrow());

  // project the data into optimized coordinate system
  data.project(gsr.axes);

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
  Bxy.row(NyDown) = Bx;

  // initialize vectors for some important derivatives
  VectorXd d2A_dx2 = VectorXd::Zero(Nx),
           dPt_dA  = VectorXd::Zero(Nx),
           d2A_dy2 = VectorXd::Zero(Nx);

  // initialize differentiator
  Differentitor differentiator;

  LOG4CPLUS_DEBUG(logger, NyUp << " steps up and " << NyDown <<
                          " steps down with step size " << dy);
  // reconstruct the upper part of the map using recursive solver
  LOG4CPLUS_DEBUG(logger, "reconstructint the upper part of the map");
  for (int i = 1; i <= NyUp; i++) {
    // 2nd derivative of A by x using Holoborodko2 filter, smoothed with
    // weighted average prior to differenting
    d2A_dx2 = differentiator.Holoborodko2(7,
      Curve::weightedAverage(Axy.row(NyDown+i-1), 1-double(abs(i)-1)/Ny/3), dx);
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
                                               1-double(abs(i))/Ny/3);
    // wtite the next row of Bx
    Bxy.row(NyDown+i) = Bxy.row(NyDown+i-1)+d2A_dy2.transpose()*dy;
    // smooth it with weighted average
    Bxy.row(NyDown+i) = Curve::weightedAverage(Bxy.row(NyDown+i),
                                               1-double(abs(i))/Ny/3);
  }

  // reconstruct the lower part of teh map using recursive solver,
  // everything is the same as for the upper part
  LOG4CPLUS_DEBUG(logger, "reconstructint the lower part of the map");
  for (int i = -1; i >= -NyDown; i--) {
    d2A_dx2 = differentiator.Holoborodko2(7,
      Curve::weightedAverage(Axy.row(NyDown+i+1), 1-double(abs(i)-1)/Ny/3), dx);
    for (int k = 0; k < Nx; k++) {
      // polyexp fitting
      dPt_dA(k) = APtFit.df(Axy(NyDown+i+1, k));
    }
    d2A_dy2 = -d2A_dx2-GSL_CONST_MKSA_VACUUM_PERMEABILITY*dPt_dA;
    Axy.row(NyDown+i) = Axy.row(NyDown+i+1)-Bxy.row(NyDown+i+1)*dy+
                        d2A_dy2.transpose()*pow(dy, 2)/2;
    Axy.row(NyDown+i) = Curve::weightedAverage(Axy.row(NyDown+i),
                                               1-double(abs(i))/Ny/3);
    Bxy.row(NyDown+i) = Bxy.row(NyDown+i+1)-d2A_dy2.transpose()*dy;
    Bxy.row(NyDown+i) = Curve::weightedAverage(Bxy.row(NyDown+i),
                                               1-double(abs(i))/Ny/3);
  }

  gsr.Axy = Axy; // save vector potential

  LOG4CPLUS_DEBUG(logger, "estimating Bz all over the map");

  // initialize the Bz(A) curve
  Curve ABzCurve(A, Bz);

  // fit it with exponent
  ExpFit ABzFit(Nx, A.data(), Bz.data());
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


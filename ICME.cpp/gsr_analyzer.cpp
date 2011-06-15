
// project headers
#include "gsr_analyzer.h"
#include "dht_analyzer.h"
#include "mva_analyzer.h"
#include "gsr_curve.h"
#include "gnuplot.h"
#include "gsl_fit_poly.h"
#include "gsl_fit_epe.h"
#include "gsl_fit_exp.h"
#include "differentiator.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_fit.h>
// standard headers
#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

// main trigger to perform GSR analysis of the event
void GsrAnalyzer::analyze(Event& event) {
  DhtAnalyzer dht; // initialize dHT analyzer
  MvaAnalyzer mva; // initialize MVA analyzer
  GsrResults  gsr; // structure for GSR results

  dht.analyze(event); // carry dHT analysis for the event

  mva.analyzePmvab(event); // carry projected MVA anaysis to get initial axes

  // make a run of axes searching algorithm, save the results in the run
  // structure
  GsrRun run = loopAxes(event, 0, 1, 90, 0, 1, 360);
  // quaternions, needed to turn to optimized axes
  Quaterniond qTheta;
  Quaterniond qPhi;

  qTheta = AngleAxisd(run.optTheta*M_PI/180, event.pmvab().axes.y);
  qPhi = AngleAxisd(run.optPhi*M_PI/180, event.pmvab().axes.z);
  run.axes.z = (qPhi*(qTheta*event.pmvab().axes.z)).normalized();
  run.axes.x = (event.dht().Vht.dot(run.axes.z)*run.axes.z-
                event.dht().Vht).normalized();
  run.axes.y = run.axes.z.cross(run.axes.x).normalized();
  run.curve = GsrCurve(event, run.axes);
  run.curve.initBranches("extremums").computeResidue();

  computeMap(event, run);
  gsr.runs.push_back(run);

  cout << gsr.runs[0].optTheta << ' ' << gsr.runs[0].optPhi << endl;

  ofstream myfile;

  myfile.open ("./rm_original.dat");
  myfile << gsr.runs[0].originalResidue << endl;
  myfile.close();

  myfile.open ("./rm_combined.dat");
  myfile << gsr.runs[0].combinedResidue << endl;
  myfile.close();
}

// run one loop of axes search algorithm
GsrRun GsrAnalyzer::loopAxes(Event& event,
                             double minTheta, double dTheta, double maxTheta,
                             double minPhi,   double dPhi,   double maxPhi)
{
  double theta, phi; // temporary theta and phi angles
  Axes axes; // temporary coordinate system
  GsrRun run; // holds the information about the current run

  // save axes limits in run structure
  run.minTheta = minTheta;
  run.dTheta   = dTheta;
  run.maxTheta = maxTheta;
  run.minPhi   = minPhi;
  run.dPhi     = dPhi;
  run.maxPhi   = maxPhi;

  // define residue matrix size
  run.originalResidue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                       int((maxPhi-minPhi)/dPhi)+1);
  run.combinedResidue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                       int((maxPhi-minPhi)/dPhi)+1);
  run.branchLength    = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                                       int((maxPhi-minPhi)/dPhi)+1);

  // temporary quaternions for making axes rotations
  Quaterniond qTheta;
  Quaterniond qPhi;

  int i, k; // angle counters

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
      run.originalResidue(i,k) = curve.originalResidue();
      run.combinedResidue(i,k) = curve.combinedResidue();
      run.branchLength(i,k) = curve.branchLength();
      phi += dPhi; // make a step in phi
      k++; // move to the next column
    } // end iteration through phi angles
    theta += dTheta; // make a step in theta step
    i++; // move to the next row
  } // end iteration through theta angles

  // searching for the minimum residue direction
  int iTheta, iPhi; // index of optimal angles
  // get indices
//  run.originalResidue.minCoeff(&iTheta, &iPhi);
  run.combinedResidue.minCoeff(&iTheta, &iPhi);
  // translate indices into optimal angles
  run.optTheta = minTheta + iTheta*dTheta;
  run.optPhi = minPhi + iPhi*dPhi;

  return run; // return
}

// compute magnetic field map for a given GSR run
GsrRun& GsrAnalyzer::computeMap(Event& event, GsrRun& run) {
  Data data(event.dataNarrow());

  data.project(run.axes);

  double dx = -event.dht().Vht.dot(run.axes.x)*event.config().samplingInterval*
               run.curve.size()/event.config().Nx;
  double dy = event.config().ratio*dx;
  int Ny = floor((event.config().maxY-event.config().minY)*
           GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy);
  dy = (event.config().maxY-event.config().minY)*
       GSL_CONST_MKSA_ASTRONOMICAL_UNIT/Ny;
  int NyUp = ceil(event.config().maxY*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy);
  int NyDown = -floor(event.config().minY*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy);
  Ny = NyDown + 1 + NyUp;
  VectorXd Y = VectorXd::Zero(Ny);
  for (int i = -NyDown; i <= NyUp; i++) {
    Y(NyDown+i) = i*dy;
  }
  int Nx = event.config().Nx;

  Curve curveIn(run.curve.branches()[0]);
  Curve curveOut(run.curve.branches()[1]);

//  Gnuplot gp1("gnuplot -persist");
//	gp1 << "p '-' w p t 'PtIn(A)', '-' w p t 'PtOut(A)'\n";
//  gp1.send(curveIn).send(curveOut);

  curveIn.resample(Nx);
  curveOut.resample(Nx);

//  Gnuplot gp2("gnuplot -persist");
//	gp2 << "p '-' w p t 'PtIn(A)', '-' w p t 'PtOut(A)'\n";
//  gp2.send(curveIn).send(curveOut);

//  Gnuplot gp3("gnuplot -persist");
//  gp3 << "p '-' w l t 'Pt(A)'\n";
//  gp3.send((Curve)run.curve);

  double xArrAll[2*Nx], yArrAll[2*Nx];

  for (int i = 0; i < 2*Nx; i++) {
    if (i < Nx) {
      xArrAll[i] = curveIn.cols().x(i);
      yArrAll[i] = curveIn.cols().y(i);
    } else {
      xArrAll[i] = curveOut.cols().x(i-Nx);
      yArrAll[i] = curveOut.cols().y(i-Nx);
    }
  }

  size_t p[2*Nx]; // permuations arrays
  gsl_sort_index(p, xArrAll, 1, 2*Nx); // sort curve data by ASC X, get permutations
  gsl_permute(p, xArrAll, 1, 2*Nx); // apply permutations
  gsl_permute(p, yArrAll, 1, 2*Nx); // apply permutations

  VectorXd xAll = VectorXd::Zero(2*Nx),
           yAll = VectorXd::Zero(2*Nx);

  int k = 0;
  for (int i = 0; i < 2*Nx; i++) {
    if (i-1 < 2*Nx && xArrAll[i] == xArrAll[i+1]) {
      continue;
    } else if (i > 0 && xArrAll[i] == xArrAll[i-1]) {
      xAll(k) = xArrAll[i];
      yAll(k) = (yArrAll[i-1]+yArrAll[i])/2;
    } else {
      xAll(k) = xArrAll[i];
      yAll(k) = yArrAll[i];
    }
    k++;
  }
  int nAll = k;
  xAll.conservativeResize(nAll);
  yAll.conservativeResize(nAll);

  Curve curveAll(xAll, yAll);

  double c0, c1, cov00, cov01, cov11, chi2;
  gsl_fit_linear(xArrAll, 1, yArrAll, 1, nAll, &c0, &c1,
                 &cov00, &cov01, &cov11, &chi2);

  int slope = (c1 > 0 ? 1 : -1);

  double pCoeff[event.config().order+1];
  gsl_fit_poly(nAll, event.config().order, xArrAll, yArrAll, pCoeff);

  VectorXd yAllPoly = VectorXd::Zero(nAll);
  for (int i = 0; i < nAll; i++) {
    yAllPoly(i) += gsl_fit_poly_eval(xAll[i], pCoeff, event.config().order);
  }

  Curve curveAllPoly(xAll, yAllPoly);

//  Gnuplot gp4("gnuplot -persist");
//  gp4 << "p '-' w p t 'Pt(A)', '-' w l t 'Pt(A) fit'\n";
//  gp4.send(curveAll).send(curveAllPoly);

  double dpCoeff[event.config().order];
  for (int i = 0; i <= event.config().order-1; i++) {
    dpCoeff[i] = pCoeff[i+1]*(i+1);
  }

  double ddpCoeff[event.config().order-1];
  for (int i = 0; i <= event.config().order-2; i++) {
    ddpCoeff[i] = dpCoeff[i+1]*(i+1);
  }

  double xb, xc;
  VectorXd xTmp;
  int xbIndex;
  if (slope > 0) {
    xb = curveAll.cols().x.minCoeff(&xbIndex);
    xc = curveAll.cols().x.maxCoeff();
    xTmp = VectorXd::LinSpaced(1000, curveAll.cols().x(0)-150,
                                     curveAll.cols().x(nAll-1)+5);
  } else {
    xb = curveAll.cols().x.maxCoeff(&xbIndex);
    xc = curveAll.cols().x.minCoeff();
    xTmp = VectorXd::LinSpaced(1000, curveAll.cols().x(0)-5,
                                     curveAll.cols().x(nAll-1)+150);
  }

  int iTmp;
  while (slope > 0 && gsl_fit_poly_eval(xb, dpCoeff,
                                        event.config().order-1) < 0 ||
         slope < 0 && gsl_fit_poly_eval(xb, dpCoeff,
                                        event.config().order-1) > 0)
  {
    if (slope > 0) {
      xb = curveAll.cols().x(++xbIndex);
    } else {
      xb = curveAll.cols().x(--xbIndex);
    }
  }

  double e2 = gsl_fit_poly_eval(xb, dpCoeff, event.config().order-1)/
              gsl_fit_poly_eval(xb, pCoeff, event.config().order);
  double e1 = gsl_fit_poly_eval(xb, pCoeff, event.config().order)/
              exp(e2*xb);
  double eCoeff[] = {e1, e2};
  double E2 = gsl_fit_poly_eval(xc, ddpCoeff, event.config().order-2)/
              gsl_fit_poly_eval(xc, dpCoeff, event.config().order-1);
  double E1 = gsl_fit_poly_eval(xc, dpCoeff, event.config().order-1)/
              E2/exp(E2*xc);
  double E3 = gsl_fit_poly_eval(xc, pCoeff, event.config().order)-
              E1*exp(E2*xc);
  double ECoeff[] = {E1, E2, E3};
  double coeff0[] = {1, -slope, 1, slope};
//  double coeff[4];
  double coeff[3];

//  gsl_fit_epe(nAll, xArrAll, yArrAll, xc, xb, pCoeff, event.config().order,
//              eCoeff, ECoeff, coeff0, coeff);

  gsl_fit_exp(nAll, xArrAll, yArrAll, coeff);

  VectorXd yTmp = VectorXd::Zero(1000);
  for (int i = 0; i < 1000; i++) {
//    yTmp(i) = gsl_fit_epe_eval_f(xTmp[i], xc, xb, coeff, pCoeff,
//                                 event.config().order, eCoeff, ECoeff);
    yTmp(i) = gsl_fit_exp_eval_f(xTmp[i], coeff);
  }

  Curve curveTmp(xTmp, yTmp);

  Gnuplot gp5("gnuplot -persist");
  gp5 <<
    "p '-' w p t 'Pt(A) in', '-' w p t 'Pt(A) out', '-' w l t 'Pt(A) fit'\n";
  gp5.send(curveIn).send(curveOut).send(curveTmp);

  VectorXd dyTmp = VectorXd::Zero(1000);
  for (int i = 0; i < 1000; i++) {
//    dyTmp(i) = gsl_fit_epe_eval_df(xTmp[i], xc, xb, coeff, pCoeff, dpCoeff,
//                                   event.config().order, eCoeff, ECoeff);
    dyTmp(i) = gsl_fit_exp_eval_df(xTmp[i], coeff);
  }

  Curve curveDerTmp(xTmp, dyTmp);

  Gnuplot gp6("gnuplot -persist");
  gp6 << "p '-' w l t 'dPt/dA'\n";
  gp6.send(curveDerTmp);

  VectorXd X = VectorXd::Zero(Nx);
  for (int i = 0; i < Nx; i++) {
    X(i) = i*dx;
  }

  VectorXd A(run.curve.cols().x);
  VectorXd Bx(data.cols().Bx);
  VectorXd By(data.cols().By);
  VectorXd Bz(data.cols().Bz);
  VectorXd Pth(data.cols().Pth);

  Curve::resample(A, Nx);
  Curve::resample(Bx, Nx);
  Curve::resample(By, Nx);
  Curve::resample(Bz, Nx);
  Curve::resample(Pth, Nx);

  MatrixXd Axy = MatrixXd::Zero(Ny, Nx),
           Bxy = MatrixXd::Zero(Ny, Nx),
           Bzy = MatrixXd::Zero(Ny, Nx);

  Axy.row(NyDown) = A;
  Bxy.row(NyDown) = Bx;

  VectorXd d2A_dx2 = VectorXd::Zero(Nx),
           dPt_dA  = VectorXd::Zero(Nx),
           d2A_dy2 = VectorXd::Zero(Nx);

  Differentitor differentiator;

  cout << NyUp << " steps up and " << NyDown << " steps down with step size " << dy << endl;

  for (int i = 1; i <= NyUp; i++) {
    d2A_dx2 = differentiator.Holoborodko2(7,
      Curve::weightedAverage(Axy.row(NyDown+i-1), 1-(abs(i)-1)/Ny/3), dx);
//    d2A_dx2 = differentiator.Holoborodko2(7, Axy.row(NyDown+i-1), dx);
    for (int k = 0; k < Nx; k++) {
      dPt_dA(k) = gsl_fit_epe_eval_df(Axy(NyDown+i-1, k), xc, xb, coeff,
        pCoeff, dpCoeff, event.config().order, eCoeff, ECoeff);
    }
    d2A_dy2 = -d2A_dx2-GSL_CONST_MKSA_VACUUM_PERMEABILITY*dPt_dA;
    Axy.row(NyDown+i) = Axy.row(NyDown+i-1)+Bxy.row(NyDown+i-1)*dy+
                        d2A_dy2.transpose()*pow(dy, 2)/2;
    Axy.row(NyDown+i) = Curve::weightedAverage(Axy.row(NyDown+i),
                                               1-abs(i)/Ny/3);
    Bxy.row(NyDown+i) = Bxy.row(NyDown+i-1)+d2A_dy2.transpose()*dy;
    Bxy.row(NyDown+i) = Curve::weightedAverage(Bxy.row(NyDown+i),
                                               1-abs(i)/Ny/3);
  }

  for (int i = -1; i >= -NyDown; i--) {
    d2A_dx2 = differentiator.Holoborodko2(7,
      Curve::weightedAverage(Axy.row(NyDown+i+1), 1-(abs(i)-1)/Ny/3), dx);
//    d2A_dx2 = differentiator.Holoborodko2(7, Axy.row(NyDown+i+1), dx);
    for (int k = 0; k < Nx; k++) {
      dPt_dA(k) = gsl_fit_epe_eval_df(Axy(NyDown+i+1, k), xc, xb, coeff,
        pCoeff, dpCoeff, event.config().order, eCoeff, ECoeff);
    }
    d2A_dy2 = -d2A_dx2-GSL_CONST_MKSA_VACUUM_PERMEABILITY*dPt_dA;
    Axy.row(NyDown+i) = Axy.row(NyDown+i+1)-Bxy.row(NyDown+i+1)*dy+
                        d2A_dy2.transpose()*pow(dy, 2)/2;
    Axy.row(NyDown+i) = Curve::weightedAverage(Axy.row(NyDown+i),
                                               1-abs(i)/Ny/3);
    Bxy.row(NyDown+i) = Bxy.row(NyDown+i+1)-d2A_dy2.transpose()*dy;
    Bxy.row(NyDown+i) = Curve::weightedAverage(Bxy.row(NyDown+i),
                                               1-abs(i)/Ny/3);
  }

//  Curve curveABz(run.curve.cols().x, data.cols().Bz);
  Curve curveABz(A, Bz);

  double e[3];
  gsl_fit_exp(Nx, A.data(), Bz.data(), e);
  VectorXd BzFit = VectorXd::Zero(Nx);
  for (int i = 0; i < Nx; i++) {
    BzFit(i) = gsl_fit_exp_eval_f(A(i), e);
  }
  Curve curveABzFit(A, BzFit);
  curveABzFit.sort().unique();

  Gnuplot gp7("gnuplot -persist");
  gp7 << "p '-' w p t 'Bz(A)', '-' w l t 'Bz(A) fitted'\n";
  gp7.send(curveABz).send(curveABzFit);

  for (int i = -NyDown; i<= NyUp; i++) {
    for (int k = 0; k < Nx; k++) {
      Bzy(NyDown+i, k) = gsl_fit_exp_eval_f(Axy(NyDown+i, k), e);
    }
  }

  ofstream myfile;
  myfile.open ("./map.dat");
  myfile << Bzy << endl;
  myfile.close();

  myfile.open ("./a.dat");
  myfile << Axy << endl;
  myfile.close();
}


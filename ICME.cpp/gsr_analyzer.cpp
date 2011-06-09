
// project headers
#include "gsr_analyzer.h"
#include "dht_analyzer.h"
#include "mva_analyzer.h"
#include "gsr_curve.h"
#include "gnuplot.h"
#include "gsl_fit_poly.h"
#include "gsr_fit.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
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
  Ny = ceil(event.config().maxY*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy)-
       floor(event.config().minY*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/dy)+1;
  VectorXd Y = VectorXd::Zero(Ny);
  for (int i = 0; i < Ny; i++) {
    Y(i) = ceil(event.config().maxY/dy)+i*dy;
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
                                     curveAll.cols().x(nAll-1)+50);
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

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;

  gsl_multifit_function_fdf f;
  double sigma[nAll];
  for (int i = 0; i < nAll; i++) sigma[i] = 1;
  gsr_fit_data params;
  params.x = xArrAll;
  params.y = yArrAll;
  params.p = pCoeff;
  params.e = eCoeff;
  params.E = ECoeff;
  params.order = event.config().order;
  params.xc = xc;
  params.xb = xb;
  params.sigma = sigma;
  params.n = nAll;

  f.f = &gsr_fit_f;
  f.df = &gsr_fit_df;
  f.fdf = &gsr_fit_fdf;
  f.n = nAll;
  f.p = 4;
  f.params = &params;

  gsl_vector *coeff0 = gsl_vector_alloc(4);
  gsl_vector_set(coeff0, 0, 1);
  gsl_vector_set(coeff0, 1, -slope);
  gsl_vector_set(coeff0, 2, 1);
  gsl_vector_set(coeff0, 3, slope);

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, nAll, 4);
  gsl_multifit_fdfsolver_set(s, &f, coeff0);

  int i = 0;
  do {
    i++;
    status = gsl_multifit_fdfsolver_iterate(s);

//    if (status) break;

    status = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
  } while (status == GSL_CONTINUE && i < 5000);

  cout << i << endl;
  cout << gsl_vector_get(s->x, 0) << endl;
  cout << gsl_vector_get(s->x, 1) << endl;
  cout << gsl_vector_get(s->x, 2) << endl;
  cout << gsl_vector_get(s->x, 3) << endl;

  double coeff[4];
  coeff[0] = gsl_vector_get(s->x, 0);
  coeff[1] = gsl_vector_get(s->x, 1);
  coeff[2] = gsl_vector_get(s->x, 2);
  coeff[3] = gsl_vector_get(s->x, 3);

  gsl_multifit_fdfsolver_free (s);

  VectorXd yTmp = VectorXd::Zero(1000);
  for (int i = 0; i < 1000; i++) {
    yTmp(i) = gsr_fit_eval_f(xTmp[i], xc, xb, coeff, pCoeff,
                             event.config().order, eCoeff, ECoeff);
  }

  Curve curveTmp(xTmp, yTmp);

  Gnuplot gp5("gnuplot -persist");
  gp5 << "p '-' w p t 'Pt(A) in', '-' w p t 'Pt(A) out', '-' w l t 'Pt(A) fit'\n";
  gp5.send(curveIn).send(curveOut).send(curveTmp);

  VectorXd dyTmp = VectorXd::Zero(1000);
  for (int i = 0; i < 1000; i++) {
    dyTmp(i) = gsr_fit_eval_df(xTmp[i], xc, xb, coeff, pCoeff, dpCoeff,
                               event.config().order, eCoeff, ECoeff);
  }

  Curve curveDerTmp(xTmp, dyTmp);

  Gnuplot gp6("gnuplot -persist");
  gp6 << "p '-' w l t 'dPt/dA'\n";
  gp6.send(curveDerTmp);

  VectorXd X = VectorXd::Zero(Nx);
  for (int i = 0; i < Nx; i++) {
    X(i) = i*dx;
  }

  MatrixXd Axy  = MatrixXd::Zero(Nx, Ny),
           Bxxy = MatrixXd::Zero(Nx, Ny),
           Bzzy = MatrixXd::Zero(Nx, Ny);
}


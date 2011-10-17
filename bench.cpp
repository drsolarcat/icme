
// project headers
#include "event.h"
#include "data.h"
#include "config.h"
#include "gsr_analyzer.h"
#include "plotter.h"
#include "my_time.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_const_mksa.h>
#include <log4cplus/logger.h>
#include <log4cplus/configurator.h>
// standard headers
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace Eigen;
using namespace log4cplus;

double sinc(double);

int main(int argc, char **argv) {

  // configure the logger
  BasicConfigurator loggerConfig;
  loggerConfig.configure();

  // get logger instance
  Logger logger = Logger::getInstance("main");

  // setting the log level
  logger.setLogLevel(ALL_LOG_LEVEL);

  if (argc != 4) {
    cout << "Not enough parameters" << endl;
    exit(0);
  } else {
    double alpha = atof(argv[1]);
    double delta = atof(argv[2]);
    double y0    = atof(argv[3]);
    double gamma = alpha/sqrt(1+pow(alpha,2));
    double zeta0 = 1/M_PI*log(1+2*gamma+2*sqrt(gamma*(gamma+1)));
    double R0    = zeta0*M_PI;
    double x0    = M_PI*cos(delta*M_PI/180);
    double y0p   = x0*tan(delta*M_PI/180)-y0*R0;

    const int nx = 300;
    VectorXd x = VectorXd::LinSpaced(nx, -R0, 1.2*R0);
    const int ny = 301;
    int iy0 = (ny-1)/2;
    VectorXd y = VectorXd::Zero(ny);
    y.head(iy0+1) = VectorXd::LinSpaced(iy0+1, -R0, 0);
    y.tail(iy0+1) = VectorXd::LinSpaced(iy0+1, 0, R0);
    Matrix2d R;
    R << cos(delta*M_PI/180),  sin(delta*M_PI/180),
         sin(delta*M_PI/180), -cos(delta*M_PI/180);
    Matrix<double, nx, ny> xx, yy, A, Bzz;
    double Amin = log(-alpha+sqrt(1+pow(alpha,2)));
    double Amax = log( alpha+sqrt(1+pow(alpha,2)));

    for (int i = 0; i < nx; i++) {
      for (int k = 0; k < ny; k++) {
        xx(i,k)  = R(0,0)*(x(i)+x0)+R(0,1)*(y(k)+y0p);
        yy(i,k)  = R(1,0)*(x(i)+x0)+R(1,1)*(y(k)+y0p);
        A(i,k)   = log(alpha*cos(xx(i,k))+cosh(yy(i,k))*sqrt(1+pow(alpha,2)));
        Bzz(i,k) = exp(-A(i,k))*pow(sinc((A(i,k)-Amin)/(Amax-Amin)),2);
      }
    }

    double BzzMin = Bzz.minCoeff();
    double BzzMax = Bzz.maxCoeff();

    Matrix<double, nx, 1> Bz = Bzz.col(iy0);
    Matrix<double, nx, 1> A0 = A.col(iy0);
    Matrix<double, nx, 1> Bx = (A.col(iy0+1)-A.col(iy0-1))/(y(iy0+1)-y(iy0-1));
    Matrix<double, nx, 1> By;
    By.segment(1, nx-2) = (A0.segment(2, nx-2)-A0.segment(0, nx-2)).array()/
                          (x.segment(2, nx-2)-x.segment(0, nx-2)).array();
    By(0)    = (A0(1)-A0(0))/(x(1)-x(0));
    By(nx-1) = (A0(nx-1)-A0(nx-2))/(x(nx-1)-x(nx-2));

    Plotter plotter(false, "./res");
    cout << A.minCoeff() << ' ' << A.maxCoeff() << endl;
    cout << Bzz.minCoeff() << ' ' << Bzz.maxCoeff() << endl;
    cout << Bx.minCoeff() << ' ' << Bx.maxCoeff() << endl;
    cout << By.minCoeff() << ' ' << By.maxCoeff() << endl;
    cout << Bz.minCoeff() << ' ' << Bz.maxCoeff() << endl;
    plotter.plotGsrMagneticMap(A.transpose(), Bzz.transpose()/1e9,
                               x/10/R0*GSL_CONST_MKSA_ASTRONOMICAL_UNIT,
                               y/10/R0*GSL_CONST_MKSA_ASTRONOMICAL_UNIT, Amax);

    double Vmc = 400*1e3; // MC speed in meters per second
    // sampling interval
    int sampling = floor((x(nx-1)-x(0))/10/R0*GSL_CONST_MKSA_ASTRONOMICAL_UNIT/Vmc/nx);

    cout << sampling << endl;

    ConfigRow configRow = {
      1, sampling, 15, 2, 'g', true, false, false, false, false, false, false,
      "BENCH", "polyexp", 0.05, -0.1, 0.1, 5, 5, 0, 90, 0, 360,
      My::Time(), My::Time()};
    DataVectors dataVectors = {
      VectorXi::Zero(nx), VectorXi::Zero(nx), VectorXi::Zero(nx),
      VectorXi::Zero(nx), VectorXi::Zero(nx), VectorXi::Zero(nx),
      sqrt(pow(Bx.array(),2)+pow(By.array(),2)+pow(Bz.array(),2)).matrix()/1e9, Bx/1e9, By/1e9, Bz/1e9, // in Tesla
      VectorXd::Ones(nx)*Vmc,
      VectorXd::Ones(nx)*Vmc, VectorXd::Zero(nx), VectorXd::Zero(nx),
      VectorXd::Zero(nx), VectorXd::Zero(nx), VectorXd::Zero(nx),
      VectorXd::Zero(nx), VectorXd::Zero(nx)};
    Data data(dataVectors);
    Event *event = new Event(configRow, data, data);
    GsrAnalyzer gsr; // GSR analysis class
    gsr.analyze(*event);




    LOG4CPLUS_DEBUG(logger, "plotting the residual maps");
    // angle arrays
    VectorXd phi, theta;

    // initialize theta array
    theta = VectorXd::LinSpaced(
      ((*event).gsr().maxTheta-(*event).gsr().minTheta)/
        (*event).gsr().dTheta+1,
      (*event).gsr().minTheta,
      (*event).gsr().maxTheta);

    // initialize phi array
    phi = VectorXd::LinSpaced(
      ((*event).gsr().maxPhi-(*event).gsr().minPhi)/
        (*event).gsr().dPhi+1,
      (*event).gsr().minPhi,
      (*event).gsr().maxPhi);

    // plot the original residual map
    plotter.plotGsrResidueMap((*event).gsr().originalResidue,
                              theta, phi,
                              (*event).gsr().optTheta,
                              (*event).gsr().optPhi,
                              (configRow.toMva ? acos(abs(event->mvab().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                              (configRow.toMva ? acos(abs(event->mvab().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL),
                              (configRow.toMva ? acos(abs(event->mvub().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                              (configRow.toMva ? acos(abs(event->mvub().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL)
                              );
    // plot the combined residual map
    plotter.plotGsrResidueMap((*event).gsr().combinedResidue,
                              theta, phi,
                              (*event).gsr().optTheta,
                              (*event).gsr().optPhi,
                              (configRow.toMva ? acos(abs(event->mvab().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                              (configRow.toMva ? acos(abs(event->mvab().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL),
                              (configRow.toMva ? acos(abs(event->mvub().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                              (configRow.toMva ? acos(abs(event->mvub().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL)
                              );

    // plot Pt(A) through matplotlib
    LOG4CPLUS_DEBUG(logger, "plotting Pt(A)");
    plotter.plotGsrAPt((*event).gsr().APtInCurve,
                       (*event).gsr().APtOutCurve,
                       (*event).gsr().APtFitCurve);

    // plot dPt/dA
    LOG4CPLUS_DEBUG(logger, "plotting dPt/dA(A)");
    plotter.plotGsrAdPt((*event).gsr().AdPtFitCurve);

    // plot Bz(A)
    LOG4CPLUS_DEBUG(logger, "plotting Bz(A)");
    plotter.plotGsrABz((*event).gsr().ABzCurve,
                       (*event).gsr().ABzFitCurve);

    LOG4CPLUS_DEBUG(logger, "plotting magnetic field map");
    // plot magnetic field map
    plotter.plotGsrMagneticMap((*event).gsr().Axy,
                               (*event).gsr().Bz,
                               (*event).gsr().X,
                               (*event).gsr().Y,
                               (*event).gsr().Ab);
  }

  return 0;
}

double sinc(double x) {
  return sin(M_PI*x)/M_PI/x;
}


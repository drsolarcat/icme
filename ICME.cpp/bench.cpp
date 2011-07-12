
// project headers
#include "event.h"
#include "data.h"
#include "plotter.h"
// library headers
#include <eigen3/Eigen/Dense>
// standard headers
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace Eigen;

double sinc(double);

int main(int argc, char **argv) {

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
        A(i,k)   = log(alpha*cos(xx(i,k))+cosh(yy(i,k)*sqrt(1+pow(alpha,2))));
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

    Plotter plotter;
    plotter.plotMagneticMap(A.transpose(), Bzz.transpose(), x, y);
  }

  return 0;
}

double sinc(double x) {
  return sin(M_PI*x)/M_PI/x;
}


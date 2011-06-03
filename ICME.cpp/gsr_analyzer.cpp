
#include "gsr_analyzer.h"
#include "dht_analyzer.h"
#include "mva_analyzer.h"
#include "gsr_curve.h"

#include <iostream>
#include <ctime>
#include <fstream>

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

// main trigger to perform GSR analysis of the event
void GsrAnalyzer::analyze(Event& event) {
  DhtAnalyzer dht; // initialize dHT analyzer
  MvaAnalyzer mva; // initialize MVA analyzer
  GsrResults  gsr;

  dht.analyze(event); // carry dHT analysis for the event

  mva.analyzePmvab(event); // carry projected MVA anaysis to get initial axes

  // make a run of axes searching algorithm, save the results in the run
  // structure
  gsr.runs.push_back(loopAxes(event, 0, 1, 90, 0, 1, 360));
//  gsr.runs.push_back(loopAxes(event,
//    gsr.runs[0].optTheta-5, 1, gsr.runs[0].optTheta+5,
//    gsr.runs[0].optPhi-5, 1, gsr.runs[0].optPhi+5));
//  gsr.runs.push_back(loopAxes(event, 0, 1, 90, 0, 1, 360));

  cout << gsr.runs[0].optTheta << ' ' << gsr.runs[0].optPhi << endl;
//  cout << gsr.runs[1].optTheta << ' ' << gsr.runs[1].optPhi << endl;
//  cout << gsr.runs[2].optTheta << ' ' << gsr.runs[2].optPhi << endl;





  // some test files
  Axes axes;
  Quaterniond qTheta, qPhi;
  qTheta = AngleAxisd(gsr.runs[0].optTheta*M_PI/180, event.pmvab().axes.y);
  axes.z = qTheta*event.pmvab().axes.z;
  qPhi = AngleAxisd(gsr.runs[0].optPhi*M_PI/180, event.pmvab().axes.z);
  axes.z = (qPhi*axes.z).normalized();
  axes.x = (event.dht().Vht.dot(axes.z)*axes.z-
            event.dht().Vht).normalized();
  axes.y = axes.z.cross(axes.x).normalized();
  GsrCurve curve(event, axes);
  curve.initBranches().computeResidue();
  cout << curve.cols().x(curve.minLeftIndex()) << ' ' <<
          curve.cols().x(curve.maxIndex()) << ' ' <<
          curve.cols().x(curve.minRightIndex()) << ' ' <<
          curve.residue() << ' ' << curve.combinedResidue() << endl;
  ofstream myfile;
  myfile.open ("./curve1.txt");
  myfile << curve.cols().x << endl;
  myfile << curve.cols().y << endl;
  myfile.close();


  qTheta = AngleAxisd(33*M_PI/180, event.pmvab().axes.y);
  axes.z = qTheta*event.pmvab().axes.z;
  qPhi = AngleAxisd(343*M_PI/180, event.pmvab().axes.z);
  axes.z = (qPhi*axes.z).normalized();
  axes.x = (event.dht().Vht.dot(axes.z)*axes.z-
            event.dht().Vht).normalized();
  axes.y = axes.z.cross(axes.x).normalized();
  curve = GsrCurve(event, axes);
  curve.initBranches().computeResidue();
  cout << curve.cols().x(curve.minLeftIndex()) << ' ' <<
          curve.cols().x(curve.maxIndex()) << ' ' <<
          curve.cols().x(curve.minRightIndex()) << ' ' <<
          curve.residue() << ' ' << curve.combinedResidue() << endl;
  myfile.open ("./curve2.txt");
  myfile << curve.cols().x << endl;
  myfile << curve.cols().y << endl;
  myfile.close();


  qTheta = AngleAxisd(35*M_PI/180, event.pmvab().axes.y);
  axes.z = qTheta*event.pmvab().axes.z;
  qPhi = AngleAxisd(318*M_PI/180, event.pmvab().axes.z);
  axes.z = (qPhi*axes.z).normalized();
  axes.x = (event.dht().Vht.dot(axes.z)*axes.z-
            event.dht().Vht).normalized();
  axes.y = axes.z.cross(axes.x).normalized();
  curve = GsrCurve(event, axes);
  curve.initBranches().computeResidue();
  cout << curve.cols().x(curve.minLeftIndex()) << ' ' <<
          curve.cols().x(curve.maxIndex()) << ' ' <<
          curve.cols().x(curve.minRightIndex()) << ' ' <<
          curve.residue() << ' ' << curve.combinedResidue() << endl;
  myfile.open ("./curve3.txt");
  myfile << curve.cols().x << endl;
  myfile << curve.cols().y << endl;
  myfile.close();


  qTheta = AngleAxisd(35*M_PI/180, event.pmvab().axes.y);
  axes.z = qTheta*event.pmvab().axes.z;
  qPhi = AngleAxisd(323*M_PI/180, event.pmvab().axes.z);
  axes.z = (qPhi*axes.z).normalized();
  axes.x = (event.dht().Vht.dot(axes.z)*axes.z-
            event.dht().Vht).normalized();
  axes.y = axes.z.cross(axes.x).normalized();
  curve = GsrCurve(event, axes);
  curve.initBranches().computeResidue();
  cout << curve.cols().x(curve.minLeftIndex()) << ' ' <<
          curve.cols().x(curve.maxIndex()) << ' ' <<
          curve.cols().x(curve.minRightIndex()) << ' ' <<
          curve.residue() << ' ' << curve.combinedResidue() << endl;
  myfile.open ("./curve4.txt");
  myfile << curve.cols().x << endl;
  myfile << curve.cols().y << endl;
  myfile.close();

  cout << "Testing" << endl;
for (int i = 180; i < 190; i++) {
  qTheta = AngleAxisd(35*M_PI/180, event.pmvab().axes.y);
  axes.z = qTheta*event.pmvab().axes.z;
  qPhi = AngleAxisd(i*M_PI/180, event.pmvab().axes.z);
  axes.z = (qPhi*axes.z).normalized();
  axes.x = (event.dht().Vht.dot(axes.z)*axes.z-
            event.dht().Vht).normalized();
  axes.y = axes.z.cross(axes.x).normalized();
  curve = GsrCurve(event, axes);
  curve.initBranches().computeResidue();
  cout << axes.x(0) << endl;
  cout << curve.residue() << endl;
}

  myfile.open ("./rml.txt");
  myfile << gsr.runs[0].residue << endl;
  myfile.close();

  myfile.open ("./l.txt");
  myfile << gsr.runs[0].length << endl;
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
  run.residue = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                               int((maxPhi-minPhi)/dPhi)+1);
  run.length  = MatrixXd::Zero(int((maxTheta-minTheta)/dTheta)+1,
                               int((maxPhi-minPhi)/dPhi)+1);

  // temporary quaternions for making axes rotations
  Quaterniond qTheta;
  Quaterniond qPhi;

  int i, k; // angle counters

//  time_t t1 = time(NULL);

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
//      cout << axes.x.norm() << ' ' << axes.y.norm() << ' ' << axes.z.norm() << endl;
//      cout << axes.x.dot(axes.y) << ' ' << axes.y.dot(axes.z) << ' ' << axes.x.dot(axes.z) << endl;
//      cout << axes.x.cross(axes.y).dot(axes.z) << endl;
      // create Pt(A) curve in new axes
//      GsrCurve* curve = new GsrCurve(event, axes);
      curve = GsrCurve(event, axes);
      // initialize branches and compute the residue
      curve.initBranches().computeResidue();
      // save residue into a matrix
      run.residue(i,k) = curve.residue();
      run.length(i,k) = curve.branchLength();
      if (theta == 35 && phi >= 180 && phi < 190) {
        cout << axes.x(0) << endl;
        cout << curve.residue() << endl;
      }
      // deallocate the curve
//      delete curve;
      phi += dPhi; // make a step in phi
      k++; // move to the next column
    } // end iteration through phi angles
    theta += dTheta; // make a step in theta step
    i++; // move to the next row
  } // end iteration through theta angles

  // searching for the minimum residue direction
  int iTheta, iPhi; // index of optimal angles
  // get indices
  run.residue.minCoeff(&iTheta, &iPhi);
  // translate indices into optimal angles
  run.optTheta = minTheta + iTheta*dTheta;
  run.optPhi = minPhi + iPhi*dPhi;

//  time_t t2 = time(NULL);
//  cout << difftime(t2, t1) << endl;

  return run; // return
}


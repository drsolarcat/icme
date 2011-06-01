
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
//  gsr.runs.push_back(loopAxes(event, 0, 5, 90, 0, 5, 360));
//  gsr.runs.push_back(loopAxes(event,
//    gsr.runs[0].optTheta-5, 1, gsr.runs[0].optTheta+5,
//    gsr.runs[0].optPhi-5, 1, gsr.runs[0].optPhi+5));
  gsr.runs.push_back(loopAxes(event, 0, 1, 90, 0, 1, 360));

  cout << gsr.runs[0].optTheta << ' ' << gsr.runs[0].optPhi << endl;
//  cout << gsr.runs[1].optTheta << ' ' << gsr.runs[1].optPhi << endl;
//  cout << gsr.runs[2].optTheta << ' ' << gsr.runs[2].optPhi << endl;





  // some test files
  Axes axes;
  Quaterniond qTheta, qPhi;
  qTheta = AngleAxisd(gsr.runs[0].optTheta*M_PI/180, event.pmvab().axes.y);
//  qTheta = AngleAxisd(8*M_PI/180, event.pmvab().axes.y);
  axes.z = qTheta*event.pmvab().axes.z;
  qPhi = AngleAxisd(gsr.runs[0].optPhi*M_PI/180, event.pmvab().axes.z);
//  qPhi = AngleAxisd(182*M_PI/180, event.pmvab().axes.z);
  axes.z = qPhi*axes.z;
  axes.x = (event.dht().Vht.dot(axes.z)*axes.z.array()-
            event.dht().Vht.array()).matrix().normalized();
  axes.y = axes.z.cross(axes.x);

  GsrCurve curve(event, axes);

  ofstream myfile;
  myfile.open ("./curve.txt");
  myfile << curve.cols().x << endl;
  myfile << curve.cols().y << endl;
  myfile.close();

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

  i = 0; // starting from 1st row
  theta = minTheta; // initialize theta
  while (theta <= maxTheta) { // begin iteration through theta angles
    // initialize theta quaternion
    qTheta = AngleAxisd(theta*M_PI/180, event.pmvab().axes.y);
    axes.z = qTheta*event.pmvab().axes.z; // rotate z axis around PMVA y axis
    phi = minPhi; // initialize phi
    k = 0; // starting from 2nd column
    while (phi <= maxPhi) { // begin iteration through phi angles
      // initialize phi quaternion
      qPhi = AngleAxisd(phi*M_PI/180, event.pmvab().axes.z);
      axes.z = qPhi*axes.z; // rotate z axis around PMVA z axis
      // initialize x axis
      axes.x = (event.dht().Vht.dot(axes.z)*axes.z.array()-
                event.dht().Vht.array()).matrix().normalized();
      // complement with y axis
      axes.y = axes.z.cross(axes.x);
      // create Pt(A) curve in new axes
      GsrCurve* curve = new GsrCurve(event, axes);
      // initialize branches and compute the residue
      (*curve).initBranches().computeResidue();
      // save residue into a matrix
      run.residue(i,k) = curve->combinedResidue();
      run.length(i,k) = curve->branchLength();
      // deallocate the curve
      delete curve;
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


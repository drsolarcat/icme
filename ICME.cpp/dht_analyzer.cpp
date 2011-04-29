
#include "dht_analyzer.h"

#include <eigen3/Eigen/Dense>
#include <gsl/gsl_statistics.h>

#include <iostream>

using namespace std;
using namespace Eigen;

// main trigger function for the dHT analysis
void DhtAnalyzer::analyze(Event& event) {

  DhtResults dht; // initialize dHT structure for storing results of analysis
  // first run of the search loop, all range of possible speeds
  // with 10km/s step
  dht.Vht = loop(event, event.dataNarrow().cols().Vx.minCoeff(), 10,
                 event.dataNarrow().cols().Vx.maxCoeff(),
                 event.dataNarrow().cols().Vy.minCoeff(), 10,
                 event.dataNarrow().cols().Vy.maxCoeff(),
                 event.dataNarrow().cols().Vz.minCoeff(), 10,
                 event.dataNarrow().cols().Vz.maxCoeff());
  // second run of the search loop 10km/s around previously found speed
  // with 1km/s step
  dht.Vht = loop(event, dht.Vht(0)-10, 1, dht.Vht(0)+10,
                        dht.Vht(1)-10, 1, dht.Vht(1)+10,
                        dht.Vht(2)-10, 1, dht.Vht(2)+10);
  // third run of the search loop 1km/s around previously found speed
  // with 0.1km/s step
  dht.Vht = loop(event, dht.Vht(0)-1, 0.1, dht.Vht(0)+1,
                        dht.Vht(1)-1, 0.1, dht.Vht(1)+1,
                        dht.Vht(2)-1, 0.1, dht.Vht(2)+1);

  // calculate Pearson's correlation coefficient between real electrical field
  // and estimated in a system moving with deHoffmann-Teller speed
  dht.cc = corr(event, dht.Vht);

  // save dHT analysis results in the event object
  event.dht(dht);
}

// this function represents a single loop in the process of the search for the
// dHT frame speed
Vector3d DhtAnalyzer::loop(Event& event,
                           double minVx, double dVx, double maxVx,
                           double minVy, double dVy, double maxVy,
                           double minVz, double dVz, double maxVz)
{
  // length of the data
  int n = event.dataNarrow().rows().size();
  Vector3d Vf, Vht; // trial frame speed and final dHT speed
  Vector3d V, B; // real-dHT speed difference and magnetic field vectors
                 // at a single timestamp
  MatrixXd E = MatrixXd::Zero(3, n); // matrix of electrical field for all
                                     // timestamps, we want to minimize it
  double D, minD = 0; // norm of the electrical field matrix and its minimum

  Vf(0) = minVx; // initialize X component of trial dHT speed
  while (Vf(0) <= maxVx) { // begin iteration through possible Vx
    Vf(1) = minVy; // initialize Y component of trial dHT speed
    while (Vf(1) <= maxVy) { // begin iteration through possible Vy
      Vf(2) = minVz; // initialize Z component of trial dHT speed
      while (Vf(2) <= maxVz) { // begin iteration through possible Vz
        // calculate electrical field
        for (int i=0; i<n; i++) {
          V(0) = event.dataNarrow().cols().Vx(i)-Vf(0);
          V(1) = event.dataNarrow().cols().Vy(i)-Vf(1);
          V(2) = event.dataNarrow().cols().Vz(i)-Vf(2);
          B(0) = event.dataNarrow().cols().Bx(i);
          B(1) = event.dataNarrow().cols().By(i);
          B(2) = event.dataNarrow().cols().Bz(i);
          E.col(i) = V.cross(B); // write electrical field vector as a column
                                 // in electrical field matrix
        }
        // calculate its norm we are minimizing
        D = (E.row(0).array().square()+
             E.row(1).array().square()+
             E.row(2).array().square()).mean();
        if (minD == 0) minD = D; // first iteration of the loop
        // save speed if new minimal norm of electrical field is found
        if (D <= minD) {
          minD = D; // new minimal norm
          Vht = Vf; // new dHT speed
        }
        Vf(2) += dVz; // make a step in Z component of the speed
      } // end iteration through possible Vz
      Vf(1) += dVy; // make a step in Y component of the speed
    } // end iteration through possible Vy
    Vf(0) += dVx; // make a step in X component of the speed
  } // end iteration through possible Vx

  return Vht; // return estimated dHT speed
}

// correlation analysis of the deHoffmann-Teller speed
double DhtAnalyzer::corr(Event& event, Vector3d Vht) {
  Vector3d V, B, E; // temporary vectors for speed and magnetic and
                    // electrical fields

  int n = event.dataNarrow().rows().size(); // length of the data
  double c1[3*n], c2[3*n]; // vectors of all components of
                           // electrical field vectors

  // begin iteration through data to calculate electrical field
  for (int i=0; i<n; i++) {
    // fill in temporary vectors
    V(0) = event.dataNarrow().row(i).Vx;
    V(1) = event.dataNarrow().row(i).Vy;
    V(2) = event.dataNarrow().row(i).Vz;
    B(0) = event.dataNarrow().row(i).Bx;
    B(1) = event.dataNarrow().row(i).By;
    B(2) = event.dataNarrow().row(i).Bz;

    // true value
    E = V.cross(B);
    c1[3*i] = E(0);
    c1[3*i+1] = E(1);
    c1[3*i+2] = E(2);
    // in dHT frame
    E = Vht.cross(B);
    c2[3*i] = E(0);
    c2[3*i+1] = E(1);
    c2[3*i+2] = E(2);
  } // end iteration through data

  // return correlation coefficient
  return gsl_stats_correlation(c1, 1, c2, 1, 3*n);
}


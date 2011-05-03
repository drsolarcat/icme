
#include "mva_analyzer.h"

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

// main trigger function to launch all the analysis
void MvaAnalyzer::analyze(Event& event) {
  analyzeMvab(event); // MVAB
  analyzeMvub(event); // MVUB
}

// MVAB analysis
void MvaAnalyzer::analyzeMvab(Event& event) {
  event.mvab(analyzeMva(event, false)); // set MVAB results in Event object
}

// MVUB analysis
void MvaAnalyzer::analyzeMvub(Event& event) {
  event.mvub(analyzeMva(event, true)); // set MVUB results in Event object
}

// PMVAB analysis
void MvaAnalyzer::analyzePmvab(Event& event) {
  event.pmvab(analyzePmva(event, false)); // set PMVAB results in Event object
}

// PMVUB analysis
void MvaAnalyzer::analyzePmvub(Event& event) {
  event.pmvub(analyzePmva(event, true)); // set PMVUB results in Event object
}

// MVA analysis
MvaResults MvaAnalyzer::analyzeMva(const Event& event, bool unit) {
  const int n = event.dataNarrow().rows().size(); // length of data vector
  MvaResults mva; // structure for MVA results

  // normalize magnetic field if unit flag is set to true
  if (unit) {
    ArrayXd norm; // dynamic norm array
    norm.resize(n); // resize norm array
    // and calculate it
    norm = (event.dataNarrow().cols().Bx.array().square()+
            event.dataNarrow().cols().By.array().square()+
            event.dataNarrow().cols().Bz.array().square()).sqrt();
    // get variance matrix using normalized magnetic field vectors
    mva.matrix = calculateVarianceMatrix(
      (event.dataNarrow().cols().Bx.array()/norm).matrix(),
      (event.dataNarrow().cols().By.array()/norm).matrix(),
      (event.dataNarrow().cols().Bz.array()/norm).matrix());
  } else {
    // get variance matrix without normalizing magnetic field vectors
    mva.matrix = calculateVarianceMatrix(
      event.dataNarrow().cols().Bx,
      event.dataNarrow().cols().By,
      event.dataNarrow().cols().Bz);
  }

  // initialize eigen solver for adjoint matrix
  SelfAdjointEigenSolver<Matrix3d> eigensolver(mva.matrix);
  // initialize and get eigen values in ascending order
  Vector3d eigenValues = eigensolver.eigenvalues();
  // initialize and get eigen vectors
  Matrix3d eigenVectors = eigensolver.eigenvectors();
  // calculate criterion of MVA analysis - ratio between intermidiate and
  // minimum eigen values
  mva.criterion = eigenValues(1)/eigenValues(0);

  // slice matrix with eigen vectors into separate vectors
  mva.axes.x = eigenVectors.col(2);
  mva.axes.y = eigenVectors.col(1);
  mva.axes.z = eigenVectors.col(0);

  // check if magnetic field is aligned with estimated axes
  if (event.dataNarrow().cols().Bx.mean()*mva.axes.y(0)+
      event.dataNarrow().cols().By.mean()*mva.axes.y(0)+
      event.dataNarrow().cols().Bz.mean()*mva.axes.y(0) < 0) {
    // if not then inverse the coordinates
    mva.axes.x = -mva.axes.x;
    mva.axes.y = -mva.axes.y;
    mva.axes.z = -mva.axes.z;
  }

  // set eigen values in MVA structure
  mva.lambdaMin = eigenValues(0);
  mva.lambdaMed = eigenValues(1);
  mva.lambdaMax = eigenValues(2);

  return mva; // return MVA results
}

// PMVA analysis
PmvaResults MvaAnalyzer::analyzePmva(const Event& event, bool unit) {
  const int n = event.dataNarrow().rows().size(); // length of data vector
  PmvaResults pmva; // structure for PMVA results

  // normalize magnetic field if unit flag is set to true
  if (unit) {
    ArrayXd norm; // dynamic norm array
    norm.resize(n); // resize norm array
    // and calculate it
    norm = (event.dataNarrow().cols().Bx.array().square()+
            event.dataNarrow().cols().By.array().square()+
            event.dataNarrow().cols().Bz.array().square()).sqrt();
    // get variance matrix using normalized magnetic field vectors
    pmva.matrix = calculateVarianceMatrix(
      (event.dataNarrow().cols().Bx.array()/norm).matrix(),
      (event.dataNarrow().cols().By.array()/norm).matrix(),
      (event.dataNarrow().cols().Bz.array()/norm).matrix());
  } else {
    // get variance matrix without normalizing magnetic field vectors
    pmva.matrix = calculateVarianceMatrix(
      event.dataNarrow().cols().Bx,
      event.dataNarrow().cols().By,
      event.dataNarrow().cols().Bz);
  }

  // initialize normalized dHT speed
  Vector3d Vht = event.dht().Vht.normalized();
  // initialize projection matrix
  pmva.pmatrix(0,0) = 1 - pow(Vht(0), 2);
  pmva.pmatrix(0,1) = -Vht(0)*Vht(1);
  pmva.pmatrix(0,2) = -Vht(0)*Vht(2);
  pmva.pmatrix(1,0) = -Vht(1)*Vht(0);
  pmva.pmatrix(1,1) = 1 - pow(Vht(1), 2);
  pmva.pmatrix(1,2) = -Vht(1)*Vht(2);
  pmva.pmatrix(2,0) = -Vht(2)*Vht(0);
  pmva.pmatrix(2,1) = -Vht(2)*Vht(1);
  pmva.pmatrix(2,2) = 1 - pow(Vht(2), 2);

  // project the variance matrix
  pmva.pmpmatrix = pmva.pmatrix*pmva.matrix*pmva.pmatrix;

  // initialize eigen solver for adjoint matrix
  SelfAdjointEigenSolver<Matrix3d> eigensolver(pmva.pmpmatrix);
  // initialize and get eigen values in ascending order
  Vector3d eigenValues = eigensolver.eigenvalues();
  // initialize and get eigen vectors
  Matrix3d eigenVectors = eigensolver.eigenvectors();

  // maximum variance vector is the Y axis
  pmva.axes.y = eigenVectors.col(2);
  // -Vht is the X axis
  pmva.axes.x = (-Vht.array()).matrix();
  // [x,y] complements the coordinates
  pmva.axes.z = (pmva.axes.x.cross(pmva.axes.y)).normalized();

  // set eigen values in PMVA structure
  pmva.lambdaMin = eigenValues(0);
  pmva.lambdaMed = eigenValues(1);
  pmva.lambdaMax = eigenValues(2);

  return pmva; // return MVA results
}

// calculate variance matrix from 3 data vectors of unknown length
Matrix3d MvaAnalyzer::calculateVarianceMatrix(const VectorXd& Bx,
                                              const VectorXd& By,
                                              const VectorXd& Bz)
{
  Matrix3d M = Matrix3d::Zero(); // initialize 3x3 matrix with zeros
  // fill the matrix with variances
  M(0,0) = (Bx.array()*Bx.array()).matrix().mean()-Bx.mean()*Bx.mean();
  M(0,1) = (Bx.array()*By.array()).matrix().mean()-Bx.mean()*By.mean();
  M(0,2) = (Bx.array()*Bz.array()).matrix().mean()-Bx.mean()*Bz.mean();
  M(1,0) = (By.array()*Bx.array()).matrix().mean()-By.mean()*Bx.mean();
  M(1,1) = (By.array()*By.array()).matrix().mean()-By.mean()*By.mean();
  M(1,2) = (By.array()*Bz.array()).matrix().mean()-By.mean()*Bz.mean();
  M(2,0) = (Bz.array()*Bx.array()).matrix().mean()-Bz.mean()*Bx.mean();
  M(2,1) = (Bz.array()*By.array()).matrix().mean()-Bz.mean()*By.mean();
  M(2,2) = (Bz.array()*Bz.array()).matrix().mean()-Bz.mean()*Bz.mean();

  return M; // return variance matrix
}


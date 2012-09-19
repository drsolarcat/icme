
// project headers
#include "config.h"
#include "data.h"
#include "event.h"
#include "mva_analyzer.h"
#include "gsr_analyzer.h"
#include "plotter.h"
// library headers
#include <eigen3/Eigen/Dense>
#include <log4cplus/logger.h>
#include <log4cplus/configurator.h>
#include <gsl/gsl_const_mksa.h>
// standard headers
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace Eigen;
using namespace My;
using namespace log4cplus;


// this is the main driver function of the whole project
int main(int argc, char* argv[]) {

  // configure the logger
  BasicConfigurator loggerConfig;
  loggerConfig.configure();

  // get logger instance
  Logger logger = Logger::getInstance("main");

  // setting the log level
  logger.setLogLevel(ALL_LOG_LEVEL);

  // parse command line arguments
  LOG4CPLUS_INFO(logger, "reading command line arguments");
  const char eventQuality = (argc > 1 ? *argv[1] : 'g'); // event quality
  LOG4CPLUS_DEBUG(logger, "events quality to process = " << eventQuality);
  // path to config file
  const string configPath = (argc > 2 ? argv[2] : "./res/config");
  LOG4CPLUS_DEBUG(logger, "path to config file = " << configPath);
  // path to directory with data files
  const string dataDir = (argc > 3 ? argv[3] : "./res");
  LOG4CPLUS_DEBUG(logger, "path to data files = " << dataDir);
  // path to directory with results
  const string resultsDir = (argc > 4 ? argv[4] : "./res");
  LOG4CPLUS_DEBUG(logger, "path to results files = " << resultsDir);

  Config config; // config object, holds config for all events
  ostringstream dataPathStream; // srting stream for data path string
  string dataPath; // path to data file
  ostringstream eventResultsDirStream; // srting stream for results directory
  string eventResultsDir, filenameResults; // path to results directory

  // initialize analyzers
  MvaAnalyzer mva; // MVA analysis class
  GsrAnalyzer gsr; // GSR analysis class


  // read the config data from the file to the config object and then
  // filter it
  LOG4CPLUS_INFO(logger, "reading and filtering the config file");
  config.readFile(configPath).filter(eventQuality);
  LOG4CPLUS_DEBUG(logger, "number of events to analyze = " <<
                          config.rows().size());

  // start iterating through the events that should be analyzed
  LOG4CPLUS_INFO(logger, "starting to iterate through events");
  for (int iEvent = 0; iEvent < config.rows().size(); iEvent++) {

    LOG4CPLUS_INFO(logger, "analysing the event " << setfill('0') <<
      config.row(iEvent).beginTime.year() << '-' <<
      setw(2) << config.row(iEvent).beginTime.month() << '-' <<
      setw(2) << config.row(iEvent).beginTime.day() << ' ' <<
      config.row(iEvent).spacecraft);
    LOG4CPLUS_DEBUG(logger, "the event length is " <<
      setiosflags(ios::fixed) << setprecision(2) <<
      double(config.row(iEvent).endTime-config.row(iEvent).beginTime)/3600 <<
      " hours");
    LOG4CPLUS_DEBUG(logger, "sampling interval = " <<
                            config.row(iEvent).samplingInterval << " seconds");
    LOG4CPLUS_DEBUG(logger, "fitting function = " <<
                            config.row(iEvent).fittingFuntion);
    LOG4CPLUS_DEBUG(logger, "polynomial order = " <<
                            config.row(iEvent).order);
    LOG4CPLUS_DEBUG(logger, "fitting parameters [boundary, center] = [" <<
                            config.row(iEvent).fittingParameterBdr << ", " <<
                            config.row(iEvent).fittingParameterCtr << ']');
    LOG4CPLUS_DEBUG(logger, "Nx = " << config.row(iEvent).Nx << " points");
    LOG4CPLUS_DEBUG(logger, "dy/dx = " << config.row(iEvent).ratio);
    LOG4CPLUS_DEBUG(logger, "minY = " << config.row(iEvent).minY << " AU");
    LOG4CPLUS_DEBUG(logger, "maxY = " << config.row(iEvent).maxY << " AU");

    LOG4CPLUS_INFO(logger, "determining the path to the data");
    // determine the path to the data file dependant on spacecraft
    // push the path to data files
    dataPathStream << dataDir << '/';
    if (config.row(iEvent).spacecraft == "WIND") {
      dataPathStream << "wind_";
      LOG4CPLUS_DEBUG(logger, "all coordinates in GSE");
    } else if (config.row(iEvent).spacecraft == "ACE") {
      dataPathStream << "ace_";
      LOG4CPLUS_DEBUG(logger, "all coordinates in GSE");
    } else if (config.row(iEvent).spacecraft == "STA") {
      dataPathStream << "stereo_a_";
      LOG4CPLUS_DEBUG(logger, "all coordinates in RTN");
    } else if (config.row(iEvent).spacecraft == "STB") {
      dataPathStream << "stereo_b_";
      LOG4CPLUS_DEBUG(logger, "all coordinates in RTN");
    } else if (config.row(iEvent).spacecraft == "CYL") {
      dataPathStream << "cylinder_";
      LOG4CPLUS_DEBUG(logger, "all coordinates in GSE");
    } else if (config.row(iEvent).spacecraft == "TOR") {
      dataPathStream << "torus_";
      LOG4CPLUS_DEBUG(logger, "all coordinates in GSE");
    } else { // throw an error - unknown spacecraft
      LOG4CPLUS_FATAL(logger, "unknown spacecraft: " <<
                              config.row(iEvent).spacecraft);
    }
    // finalize by adding resolution
    dataPathStream << config.row(iEvent).samplingInterval << ".dat";
    // save the data path
    dataPath = dataPathStream.str();
    LOG4CPLUS_DEBUG(logger, "path to the data file = " << dataPath);
    // clear the stream
    dataPathStream.clear();
    dataPathStream.str("");

    // create Time objects for wider limits of the event, initially copy
    // existing time limits
    Time* preBeginTime = new Time(config.row(iEvent).beginTime);
    Time* postEndTime = new Time(config.row(iEvent).endTime);
    // widen time limits
    preBeginTime->add(-1, "day");
    postEndTime->add(1, "day");
    // create Data object for wider time limits
    Data* dataWide = new Data(); // create dynamic object for data
    LOG4CPLUS_INFO(logger, "processing the data");
    LOG4CPLUS_DEBUG(logger, "reading the data, filtered by a wider time period");
    dataWide->readFile(dataPath, *preBeginTime, *postEndTime);
    // create Data object for original time limits
    LOG4CPLUS_DEBUG(logger, "wide data size = " << dataWide->rows().size());
    Data* dataNarrow = new Data(*dataWide);
    // filter narrow Data object out of wider one
    LOG4CPLUS_DEBUG(logger, "filtering the data to a narrower time period");
    dataNarrow->filter(config.row(iEvent).beginTime,
                       config.row(iEvent).endTime);
    LOG4CPLUS_DEBUG(logger, "narrow data size = " << dataNarrow->rows().size());
    // create dynamic object to store all event data and results of analysis
    Event* event = new Event(config.row(iEvent), *dataWide, *dataNarrow,
                             dataDir);

    LOG4CPLUS_DEBUG(logger, "average MC speed = " << setiosflags(ios::fixed) <<
      setprecision(2) << event->dataNarrow().cols().Vp.mean()/1e3 << "km/s");
    LOG4CPLUS_DEBUG(logger, "average SW speed = " << setiosflags(ios::fixed) <<
      setprecision(2) << (event->dataWide().cols().Vp.sum()-
        event->dataNarrow().cols().Vp.sum())/
        (event->dataWide().cols().Vp.size()-
        event->dataNarrow().cols().Vp.size())/1e3 << "km/s");

    // estimate the initiation time
    Time cmeTime = (*event).config().beginTime;
    cmeTime.add(-GSL_CONST_MKSA_ASTRONOMICAL_UNIT/
                 abs((*event).dataNarrow().cols().Vx.mean()), "second");
    LOG4CPLUS_INFO(logger, "Estimated CME time: " << setfill('0') <<
      cmeTime.year() << '-' <<
      setw(2) << cmeTime.month() << '-' <<
      setw(2) << cmeTime.day() << ' ' <<
      setw(2) << cmeTime.hour() << ':' <<
      setw(2) << cmeTime.minute() << ':' <<
      setw(2) << cmeTime.second());

    // initialize the directory for results
    eventResultsDirStream << resultsDir << "/" << setfill('0') <<
      config.row(iEvent).beginTime.year() << "-" <<
      setw(2) << config.row(iEvent).beginTime.month() << "-" <<
      setw(2) << config.row(iEvent).beginTime.day() << "_" <<
      config.row(iEvent).spacecraft;
    // save the results directory
    eventResultsDir = eventResultsDirStream.str();
    LOG4CPLUS_DEBUG(logger, "results directory = " << eventResultsDir);
    // the filename for numerical results
    eventResultsDirStream << "/info.txt";
    filenameResults = eventResultsDirStream.str();
    // clear the stream
    eventResultsDirStream.clear();
    eventResultsDirStream.str("");

    // plot residue maps through Matlab
    Plotter plotter(config.row(iEvent).toSave, eventResultsDir);

    // perform MVA analysis if required
    if (config.row(iEvent).toMva) {
      LOG4CPLUS_INFO(logger, "doing MVA analysis");
      mva.analyze(*event);

      LOG4CPLUS_INFO(logger, "plotting B rotation for MVAB");
      Data dataMvab(event->dataNarrow());
      dataMvab.project((*event).mvab().axes);
      plotter.plotMvaBrot(dataMvab.cols().Bx, dataMvab.cols().By);

      LOG4CPLUS_INFO(logger, "plotting B rotation for MVUB");
      Data dataMvub(event->dataNarrow());
      dataMvub.project((*event).mvub().axes);
      plotter.plotMvaBrot(dataMvub.cols().Bx, dataMvub.cols().By);
    }

    // perform GSR analysis if required
    if (config.row(iEvent).toGsr) {

      LOG4CPLUS_INFO(logger, "doing GSR analysis");
      gsr.analyze(*event); // do GSR analysis

      LOG4CPLUS_INFO(logger, "plotting results of GSR analysis");

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
                                (config.row(iEvent).toMva ? acos(abs(event->mvab().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                                (config.row(iEvent).toMva ? acos(abs(event->mvab().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL),
                                (config.row(iEvent).toMva ? acos(abs(event->mvub().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                                (config.row(iEvent).toMva ? acos(abs(event->mvub().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL)
                                );
      // plot the combined residual map
      plotter.plotGsrResidueMap((*event).gsr().combinedResidue,
                                theta, phi,
                                (*event).gsr().optTheta,
                                (*event).gsr().optPhi,
                                (config.row(iEvent).toMva ? acos(abs(event->mvab().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                                (config.row(iEvent).toMva ? acos(abs(event->mvab().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL),
                                (config.row(iEvent).toMva ? acos(abs(event->mvub().axes.y.dot(event->pmvab().axes.z)))*180/M_PI : NULL),
                                (config.row(iEvent).toMva ? acos(abs(event->mvub().axes.y.cross(event->pmvab().axes.z).dot(event->pmvab().axes.x)))*180/M_PI+90 : NULL)
                                );

      // plot full Pt(A) through matplotlib
      LOG4CPLUS_DEBUG(logger, "plotting full Pt(A)");
      plotter.plotGsrAPtFull((*event).gsr().curve);

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
                                 (*event).gsr().Ab,
                                 (*event).gsr().Aa,
                                 (*event).gsr().Bx,
                                 (*event).gsr().By,
                                 (*event).gsr().axes);
    }

    // plot the in-situ data
    plotter.plotData(*event);

    if (config.row(iEvent).toSave) {
      ofstream infoFile(filenameResults.c_str());

      infoFile << "Event " << setfill('0') <<
        config.row(iEvent).beginTime.year() << '-' <<
        setw(2) << config.row(iEvent).beginTime.month() << '-' <<
        setw(2) << config.row(iEvent).beginTime.day() << ' ' <<
        config.row(iEvent).spacecraft << endl;

      infoFile << "The event length (red boundaries) is " <<
        setiosflags(ios::fixed) << setprecision(2) <<
        double(config.row(iEvent).endTime-config.row(iEvent).beginTime)/3600 <<
        " hours" << endl;

      infoFile << "Estimated CME time: " << setfill('0') <<
        cmeTime.year() << '-' <<
        setw(2) << cmeTime.month() << '-' <<
        setw(2) << cmeTime.day() << ' ' <<
        setw(2) << cmeTime.hour() << ':' <<
        setw(2) << cmeTime.minute() << ':' <<
        setw(2) << cmeTime.second() << endl;

      if (config.row(iEvent).toGsr) {
        infoFile << "Estimated de Hoffmann-Teller speed = [" <<
          setiosflags(ios::fixed) << setprecision(1) <<
          (*event).dht().Vht(0)/1e3 << ", " <<
          (*event).dht().Vht(1)/1e3 << ", " <<
          (*event).dht().Vht(2)/1e3 << "] km/s" << endl;

        infoFile << "Magnetic cloud axes: " <<
          setiosflags(ios::fixed) << setprecision(3) << "x[" <<
          (*event).gsr().axes.x(0) << ", " <<
          (*event).gsr().axes.x(1) << ", " <<
          (*event).gsr().axes.x(2) << "], y[" <<
          (*event).gsr().axes.y(0) << ", " <<
          (*event).gsr().axes.y(1) << ", " <<
          (*event).gsr().axes.y(2) << "], z[" <<
          (*event).gsr().axes.z(0) << ", " <<
          (*event).gsr().axes.z(1) << ", " <<
          (*event).gsr().axes.z(2) << "]" << endl;

        infoFile << "HEEQ magnetic cloud axes: " <<
          setiosflags(ios::fixed) << setprecision(3) << "x[" <<
          (*event).gsr().heeqAxes.x(0) << ", " <<
          (*event).gsr().heeqAxes.x(1) << ", " <<
          (*event).gsr().heeqAxes.x(2) << "], y[" <<
          (*event).gsr().heeqAxes.y(0) << ", " <<
          (*event).gsr().heeqAxes.y(1) << ", " <<
          (*event).gsr().heeqAxes.y(2) << "], z[" <<
          (*event).gsr().heeqAxes.z(0) << ", " <<
          (*event).gsr().heeqAxes.z(1) << ", " <<
          (*event).gsr().heeqAxes.z(2) << "]" << endl;

        infoFile << "Invariant axis in Stonyhurst coordinates " <<
          "[theta, phi] = [" << setiosflags(ios::fixed) << setprecision(1) <<
          (*event).gsr().stonyhurstTheta << ", " <<
          (*event).gsr().stonyhurstPhi << "]" << endl;

        infoFile << "GSR boundaries: " << setfill('0') <<
          (*event).gsr().beginTime.year() << '-' <<
          setw(2) << (*event).gsr().beginTime.month() << '-' <<
          setw(2) << (*event).gsr().beginTime.day() << ' ' <<
          setw(2) << (*event).gsr().beginTime.hour() << ':' <<
          setw(2) << (*event).gsr().beginTime.minute() << ':' <<
          setw(2) << (*event).gsr().beginTime.second() << " - " <<
          (*event).gsr().endTime.year() << '-' <<
          setw(2) << (*event).gsr().endTime.month() << '-' <<
          setw(2) << (*event).gsr().endTime.day() << ' ' <<
          setw(2) << (*event).gsr().endTime.hour() << ':' <<
          setw(2) << (*event).gsr().endTime.minute() << ':' <<
          setw(2) << (*event).gsr().endTime.second() << endl;

        infoFile << "Impact parameter = " << setprecision(8) <<
          (*event).gsr().ip/GSL_CONST_MKSA_ASTRONOMICAL_UNIT << "AU" << endl;

        infoFile << "Boundary value of A = " <<
          (*event).gsr().Ab << endl;
      }
    }

    LOG4CPLUS_INFO(logger, "everything is done");

    // delete dynamic objects
    delete dataWide;
    delete dataNarrow;
    delete preBeginTime;
    delete postEndTime;
    delete event;
  } // end of iteration through the events
}



// project headers
#include "config.h"
#include "data.h"
#include "event.h"
#include "mva_analyzer.h"
#include "gsr_analyzer.h"
#include "gnuplot.h"
#include "plotter.h"
// library headers
#include "engine.h"
#include <eigen3/Eigen/Dense>
#include <log4cplus/logger.h>
#include <log4cplus/configurator.h>
// standard headers
#include <string>
#include <iostream>
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
  string eventResultsDir; // path to results directory

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
    preBeginTime->add(-3, "hour");
    postEndTime->add(3, "hour");
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
    Event* event = new Event(config.row(iEvent), *dataWide, *dataNarrow);

    // initialize the directory for results
    eventResultsDirStream << resultsDir << "/" << setfill('0') <<
      config.row(iEvent).beginTime.year() << "-" <<
      setw(2) << config.row(iEvent).beginTime.month() << "-" <<
      setw(2) << config.row(iEvent).beginTime.day() << "_" <<
      config.row(iEvent).spacecraft;
    // save the results directory
    eventResultsDir = eventResultsDirStream.str();
    LOG4CPLUS_DEBUG(logger, "results directory = " << eventResultsDir);
    // clear the stream
    eventResultsDirStream.clear();
    eventResultsDirStream.str("");

    // perform GSR analysis if required
    if (config.row(iEvent).toGsr) {

      LOG4CPLUS_INFO(logger, "doing GSR analysis");
      gsr.analyze(*event); // do GSR analysis

      LOG4CPLUS_INFO(logger, "plotting results of GSR analysis");
      // plot residue maps through Matlab
      Plotter plotter(config.row(iEvent).toSave, eventResultsDir);

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
//      plotter.plotResidueMap((*event).gsr().originalResidue,
//                             theta, phi,
//                             (*event).gsr().optTheta,
//                             (*event).gsr().optPhi);

      // plot the combined residual map
      plotter.plotResidueMap((*event).gsr().combinedResidue,
                             theta, phi,
                             (*event).gsr().optTheta,
                             (*event).gsr().optPhi);

      // plot Pt(A) through matplotlib
      LOG4CPLUS_DEBUG(logger, "plotting Pt(A)");
      plotter.plotGsrAPt((*event).gsr().APtInCurve,
                         (*event).gsr().APtOutCurve,
                         (*event).gsr().APtFitCurve);

      // plot dPt/dA through Gnuplot
      LOG4CPLUS_DEBUG(logger, "plotting dPt/dA(A)");
      plotter.plotGsrAdPt((*event).gsr().AdPtFitCurve);

      // plot Bz(A) through Gnuplot
      LOG4CPLUS_DEBUG(logger, "plotting Bz(A)");
      plotter.plotGsrABz((*event).gsr().ABzCurve,
                         (*event).gsr().ABzFitCurve);

      LOG4CPLUS_DEBUG(logger, "plotting magnetic field map");
      // plot magnetic field map through Matlab
      plotter.plotMagneticMap((*event).gsr().Axy,
                              (*event).gsr().Bz,
                              (*event).gsr().X,
                              (*event).gsr().Y,
                              (*event).gsr().Ab);
    }

    // perform MVA analysis if required
    if (config.row(iEvent).toMva) {
      LOG4CPLUS_INFO(logger, "doing MVA analysis");
      mva.analyze(*event);
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


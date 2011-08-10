
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

  // set logger level to info
//  logger.setLogLevel(INFO_LOG_LEVEL);
  logger.setLogLevel(DEBUG_LOG_LEVEL);

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

  Config config; // config object, holds config for all events
  ostringstream dataPathStream; // srting stream for data path string
  string dataPath; // path to data file

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

    LOG4CPLUS_INFO(logger, "determining the path to the data");
    // determine the path to the data file dependant on spacecraft
    // push the path to data files
    dataPathStream << dataDir << '/';
    if (config.row(iEvent).spacecraft == "WIND") {
      dataPathStream << "wind_";
    } else if (config.row(iEvent).spacecraft == "ACE") {
      dataPathStream << "ace_";
    } else if (config.row(iEvent).spacecraft == "STA") {
      dataPathStream << "stereo_a_";
    } else if (config.row(iEvent).spacecraft == "STB") {
      dataPathStream << "stereo_b_";
    } else { // throw an error - unknown spacecraft
      cout << "Unknown spacecraft" << endl;
    }
    // finalize by adding resolution
    dataPathStream << config.row(iEvent).samplingInterval << ".dat";
    // save the data path
    dataPath = dataPathStream.str();
    LOG4CPLUS_DEBUG(logger, "path to the data file = " << dataPath);
    // clear the stream
    dataPathStream.clear();

    // create Time objects for wider limits of the event, initially copy
    // existing time limits
    Time* preBeginTime = new Time(config.row(iEvent).beginTime);
    Time* postEndTime = new Time(config.row(iEvent).endTime);
    // widen time limits
    preBeginTime->add(-3, "hour");
    postEndTime->add(3, "hour");
    // create Data object for wider time limits
    Data* dataWide = new Data(); // create dynamic object for data
    LOG4CPLUS_INFO(logger, "reading the data, filtered by a wider time period");
    dataWide->readFile(dataPath, *preBeginTime, *postEndTime);
    // create Data object for original time limits
    LOG4CPLUS_DEBUG(logger, "wide data size = " << dataWide->rows().size());
    Data* dataNarrow = new Data(*dataWide);
    // filter narrow Data object out of wider one
    LOG4CPLUS_INFO(logger, "filtering the data to a narrower time period");
    dataNarrow->filter(config.row(iEvent).beginTime,
                       config.row(iEvent).endTime);
    LOG4CPLUS_DEBUG(logger, "narrow data size = " << dataNarrow->rows().size());
    // create dynamic object to store all event data and results of analysis
    Event* event = new Event(config.row(iEvent), *dataWide, *dataNarrow);

    // perform GSR analysis if required
    if (config.row(iEvent).toGsr) {

      LOG4CPLUS_INFO(logger, "starting GSR analysis");
      gsr.analyze(*event); // do GSR analysis

      if (config.row(iEvent).toPlot) { // plot

        LOG4CPLUS_INFO(logger, "plotting results of GSR analysis");
        // plot residue maps through Matlab
        Plotter plotter;

        cout << "plotting residual maps... ";
        // angle arrays
        VectorXd phi, theta;

        theta = VectorXd::LinSpaced(
          ((*event).gsr().maxTheta-(*event).gsr().minTheta)/
            (*event).gsr().dTheta+1,
          (*event).gsr().minTheta,
          (*event).gsr().maxTheta);

        phi = VectorXd::LinSpaced(
          ((*event).gsr().maxPhi-(*event).gsr().minPhi)/
            (*event).gsr().dPhi+1,
          (*event).gsr().minPhi,
          (*event).gsr().maxPhi);

        plotter.plotResidueMap((*event).gsr().originalResidue,
                               theta, phi,
                               (*event).gsr().optTheta,
                               (*event).gsr().optPhi);

        plotter.plotResidueMap((*event).gsr().combinedResidue,
                               theta, phi,
                               (*event).gsr().optTheta,
                               (*event).gsr().optPhi);
        cout << "done" << endl;

        cout << "plotting Pt(A), dPt/dA(A) and Bz(A) curves... ";
        // plot Pt(A) through Gnuplot
        Gnuplot APt("gnuplot -persist");
        APt << "p '-' w p t 'Pt(A) in', "
              << "'-' w p t 'Pt(A) out', "
              << "'-' w l t 'Pt(A) fit'\n";
        APt.send((*event).gsr().APtInCurve).
            send((*event).gsr().APtOutCurve).
            send((*event).gsr().APtFitCurve);

        // plot dPt/dA through Gnuplot
        Gnuplot AdPt("gnuplot -persist");
        AdPt << "p '-' w l t 'dPt/dA fit'\n";
        AdPt.send((*event).gsr().AdPtFitCurve);

        // plot Bz(A) through Gnuplot
        Gnuplot ABz("gnuplot -persist");
        ABz << "p '-' w p t 'Bz(A)', '-' w l t 'Bz(A) fit'\n";
        ABz.send((*event).gsr().ABzCurve).
            send((*event).gsr().ABzFitCurve);
        cout << "done" << endl;

        cout << "plotting magnetic field map... ";
        // plot magnetic field map through Matlab
        plotter.plotMagneticMap((*event).gsr().Axy,
                                (*event).gsr().Bz,
                                (*event).gsr().X,
                                (*event).gsr().Y);
        cout << "done" << endl;
        cout << "done GSR plotting" << endl;
      }
    }

    // perform MVA analysis if required
    if (config.row(iEvent).toMva) {
      cout << "starting MVA analysis" << endl;
      mva.analyze(*event);
      cout << "done MVA analysis" << endl;
    }

    // delete dynamic objects
    delete dataWide;
    delete dataNarrow;
    delete preBeginTime;
    delete postEndTime;
    delete event;
  } // end of iteration through the events
}


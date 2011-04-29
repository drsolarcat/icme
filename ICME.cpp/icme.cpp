
#include "config.h"
#include "data.h"
#include "event.h"
#include "mva_analyzer.h"
#include "gsr_analyzer.h"

#include <string>
#include <iostream>

using namespace std;

// this is the main driver function of the whole project
int main() {
  const string configPath = "./res/config"; // path to config file
  // the quality of the events to be analyzed
  // TODO: it should be in command line arguments
  const char eventQuality = 'g';
  Config config; // config object, holds config for all events
  string dataPath; // path to data file, will be changed for each event
  // initialize analyzers
  MvaAnalyzer mva; // MVA analysis class
  GsrAnalyzer gsr; // GSR analysis class

  // read the config data from the file to the config object and then
  // filter it
  config.readFile(configPath).filter(eventQuality);

  // start iterating through the events that should be analyzed
  for (int iEvent = 0; iEvent < config.rows().size(); iEvent++) {

    // determine the path to the data file dependant on spacecraft
    if (config.row(iEvent).spacecraft == "WIND") {
      dataPath = "./res/wind_240.dat";
    } else if (config.row(iEvent).spacecraft == "ACE") {
      dataPath = "./res/ace_240.dat";
    } else if (config.row(iEvent).spacecraft == "STA") {
      dataPath = "./res/stereo_a_240.dat";
    } else if (config.row(iEvent).spacecraft == "STB") {
      dataPath = "./res/stereo_b_240.dat";
    } else { // throw an error - unknown spacecraft
      cout << "Unknown spacecraft" << endl;
    }

    // create Time objects for wider limits of the event, initially copy
    // existing time limits
    Time* preBeginTime = new Time(config.row(iEvent).beginTime);
    Time* postEndTime = new Time(config.row(iEvent).endTime);
    // widen time limits
    preBeginTime->add(-3, "hour");
    postEndTime->add(3, "hour");
    // create Data object for wider time limits
    Data* dataWide = new Data(); // create dynamic object for data
    dataWide->readFile(dataPath, *preBeginTime, *postEndTime);
    // create Data object for original time limits
    Data* dataNarrow = new Data(*dataWide);
    // filter narrow Data object out of wider one
    dataNarrow->filter(config.row(iEvent).beginTime,
      config.row(iEvent).endTime);

    // create dynamic object to store all event data and results of analysis
    Event* event = new Event(config.row(iEvent), *dataWide, *dataNarrow);

    // perform GSR analysis if required
    if (config.row(iEvent).toGsr) gsr.analyze(*event);

    // perform MVA analysis if required
    if (config.row(iEvent).toMva) mva.analyze(*event);

    // delete dynamic objects
    delete dataWide;
    delete dataNarrow;
    delete preBeginTime;
    delete postEndTime;
    delete event;
  } // end of iteration through the events
}


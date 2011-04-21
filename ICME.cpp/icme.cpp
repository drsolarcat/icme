
#include "config.h"

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

  // read the config data from the file to the config object and then filter it
  config.readFile(configPath).filter(eventQuality);

  // start iterating through the events that should be analyzed
  for (int iEvent = 0; iEvent < config.rows().size(); iEvent++) {
    // NEXT: get data for the analyzed event
    cout << config.row(iEvent).beginDate << endl;
  } // end of iteration through the events
}


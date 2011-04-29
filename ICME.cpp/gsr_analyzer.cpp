
#include "gsr_analyzer.h"

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

// main trigger to perform GSR analysis of the event
void GsrAnalyzer::analyze(Event& event) {
  DhtAnalyzer dht; // initialize dHT analyzer

  dht.analyze(event); // carry dHT analysis for the event
//  event.gsr(analyzeGsr(event)); // carry GSR analysis for the event
}


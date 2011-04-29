
#ifndef GSR_ANALYZER_H
#define GSR_ANALYZER_H

#include "event.h"
#include "dht_analyzer.h"

// this class is used to carry GSR analysis
class GsrAnalyzer {
  public:
    void analyze(Event&); // main trigger to carry the analysis,
                          // the results are stored in the Event object
//  private:
//    GsrResults analyzeGsr(Event&); // the main part of the analysis
};

#endif


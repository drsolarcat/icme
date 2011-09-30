
// project headers
#include "event.h"

using namespace std;

// constructor for the Event object, configuration and data are passed in
Event::Event(ConfigRow config, Data dataWide, Data dataNarrow) {
  _config = config;
  _dataWide = dataWide;
  _dataNarrow = dataNarrow;
}


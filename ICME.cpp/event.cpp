
#include "event.h"

using namespace std;

// constructor for the Event object, configuration and data are passed in
Event::Event(ConfigRow config, Data dataWide, Data dataNarrow) {
  e_config = config;
  e_dataWide = dataWide;
  e_dataNarrow = dataNarrow;
}


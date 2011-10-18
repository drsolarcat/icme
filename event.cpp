
// project headers
#include "event.h"
// standard headers
#include <string>

using namespace std;

// constructor for the Event object, configuration and data are passed in
Event::Event(ConfigRow config, Data dataWide, Data dataNarrow, string dataDir) :
  _config(config), _dataWide(dataWide), _dataNarrow(dataNarrow),
  _dataDir(dataDir)
{}



#include "time.h"

#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main() {
  int year = 1970;
  int month = 1;
  int day = 1;
  int hour = 12;
  int minute = 0;
  int second = 0;

  Time mytime(year, month, day, hour, minute, second);

  cout << fixed << setprecision(4) << mytime.matlabtime() << endl;
}


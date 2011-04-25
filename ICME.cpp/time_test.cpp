
#include "time.h"

#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main() {
  int year = 1970;
  int month = 1;
  int day = 2;
  int hour = 12;
  int minute = 0;
  int second = 0;
  int doy = 2;
  int unixtime = 129600;
  double matlabtime = 719530.5;
  const string stringtime = "1970-01-02 12:00:00";

  Time mytime1(year, month, day, hour, minute, second);
  cout << fixed << setprecision(4) << mytime1.matlabtime() << endl;

  Time mytime2(year, doy, hour, minute, second);
  cout << fixed << setprecision(4) << mytime2.matlabtime() << endl;

  Time mytime3(unixtime, "unix");
  cout << fixed << setprecision(4) << mytime3.matlabtime() << endl;

  Time mytime4(matlabtime, "matlab");
  cout << fixed << setprecision(4) << mytime4.matlabtime() << endl;

  Time mytime5(stringtime);
  cout << fixed << setprecision(4) << mytime5.matlabtime() << endl;

  Time mytime6;
  cout << fixed << setprecision(4) << mytime6.matlabtime() << endl;
}


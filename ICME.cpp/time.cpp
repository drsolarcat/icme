
#include "time.h"

#include <string>
#include <ctime>
#include <iostream>
#include <iomanip> // for setprecision() function

using namespace std;

// constructors

// construct Time object using current time
Time::Time() {
  // initialize local-UTC time shift
  initLocalUtcShift();

  // initialize ctime structure
  initTm();

  c_unixtime = time(NULL); // current Unix time
  c_tm = *gmtime(&c_unixtime); // ctime structure of the current time

  // setting Time object values for the current time
  t_year = c_tm.tm_year+1900; // exact value of year
  t_month = c_tm.tm_mon+1; // in 1..12 range
  t_day = c_tm.tm_mday; // in 1..31 range
  t_hour = c_tm.tm_hour; // in 0..23 range
  t_minute = c_tm.tm_min; // in 0..59 range
  t_second = c_tm.tm_sec; // in 0..61 range, leap seconds included
  // dropping off leap seconds
  if (t_second > 59) t_second = 59; // now in 0..59 range
  t_doy = c_tm.tm_yday+1; // we store days of year in 1..366 range

  t_unixtime = int(c_unixtime); // cast Unix time from time_t to int and save

  // recalculate ctime structure to get day of year
  c_tm = *gmtime(const_cast<time_t*>(&c_unixtime));
  // save day of year to the Time object
  t_doy = c_tm.tm_yday+1; // we want it in 1..366 range

  // calculate Matlab time
  calculateMatlabTime();
}

// construct Time object using year, month, day, hour and second
Time::Time(int year, int month, int day, int hour, int minute, int second) {
  // initialize local-UTC time shift
  initLocalUtcShift();

  // initialize ctime structure
  initTm();

  // save provided time data
  t_year = year;
  t_month = month;
  t_day = day;
  t_hour = hour;
  t_minute = minute;
  t_second = second;

  // fill it with time data we already know
  c_tm.tm_year = t_year-1900; // number of years since 1900
  c_tm.tm_mon = t_month-1; // in range 0..11
  c_tm.tm_mday = t_day;
  c_tm.tm_hour = t_hour+t_local_utc_shift; // taking timezone into account
  c_tm.tm_min = t_minute;
  c_tm.tm_sec = t_second;
  // calculate Unix time base on ctime structure
  c_unixtime = mktime(&c_tm);
  // cast it into integer and save in the Time object
  t_unixtime = int(c_unixtime);

  // recalculate ctime structure to get day of year
  c_tm = *gmtime(const_cast<time_t*>(&c_unixtime));
  // save day of year to the Time object
  t_doy = c_tm.tm_yday+1; // we want it in 1..366 range

  // calculate Matlab time
  calculateMatlabTime();
}

// construct Time object using year, day of year, hour, minute and second
Time::Time(int year, int doy, int hour, int minute, int second) {
  // initialize local-UTC time shift
  initLocalUtcShift();

  // initialize ctime strucure
  initTm();

  // save provided time data
  t_year = year;
  t_doy = doy;
  t_hour = hour;
  t_minute = minute;
  t_second = second;

  // fill ctime strucure with available time data
  c_tm.tm_year = t_year;
  c_tm.tm_yday = t_doy;
  c_tm.tm_hour = t_hour+t_local_utc_shift; // take into account the timezone
  c_tm.tm_min = t_minute;
  c_tm.tm_sec = t_second;
  // calculate Unix time
  c_unixtime = mktime(&c_tm);
  // cast it into integer
  t_unixtime = int(c_unixtime);

  // recalculate ctime structure to get day of month
  c_tm = *gmtime(const_cast<time_t*>(&c_unixtime));
  // save day of year to the Time object
  t_day = c_tm.tm_mday;

  // calculate Matlab time
  calculateMatlabTime();
}

// construct Time object using formatted string representation of time
Time::Time(string dateTimeString) {}

// construct Time object using Unix, Matlab or MySQL timestamp
Time::Time(int timestamp, string type) {
  // check the type of the timestamp provided
  if (type == "unix") { // if Unix timestamp is provided

  } else if (type == "matlab") { // if Matlab timestamp is provided

  } else if (type == "mysql") { // if MySQL timestamp is provided

  } else { // throw an error - unknown type of the timestamp

  }
}

// other methods

// add "amount" of "type" (year, month, day, hour, minute, second) time to the
// Time object
Time& Time::add(int amount, string type) {}

// initialize shift between local and UTC time
void Time::initLocalUtcShift() {
  // current time
  time_t c_time = time(NULL);
  // initialize ctyme structures for local and UTC times
  tm * c_tm_local, * c_tm_utc;
  // evaluate local time structure
  c_tm_local = localtime(&c_time);
  // get Unix time for the local time
  int c_sec_local = mktime(c_tm_local);
  // evaluate UTC time structure
  c_tm_utc = gmtime(&c_time);
  // get Unix time for the UTC time
  int c_sec_utc = mktime(c_tm_utc);
  // calculate and save time shift between local and UTC time in hours
  t_local_utc_shift = (c_sec_local-c_sec_utc)/3600;
}

// initialize ctime strucure with zeros
void Time::initTm() {
  tm c_tm = {0};
}

// calculate Matlab timestamp from ctime strucure and Unix timestamp
void Time::calculateMatlabTime() {
  // calculating the Matlab timestamp
  // number of seconds between absolute zero time (0000-01-01 00:00:00) is
  // too large, so we will calculate this difference for Unix zero time
  // create ctime structure for Unix epoch zero time
  tm c_tm0 = {0};
  c_tm0.tm_year = 1970-1900; // number of years since 1900
  c_tm0.tm_mon = 0;
  c_tm0.tm_mday = 1;
  c_tm0.tm_hour = 0+t_local_utc_shift; // taking timezone into account
  c_tm0.tm_min = 0;
  c_tm0.tm_sec = 0;
  // get Unix time for zero time
  time_t c_unixtime0 = mktime(&c_tm0);
  // calculate elapsed seconds since zero time
  double c_sec0 = difftime(c_unixtime, c_unixtime0);
  // calculate elapsed days since zero time
  // 719529 is the number of days from 0000-01-01 to 1970-01-01
  // 86400 is the number of seconds in one day
  double c_matlabtime = double(719529)+c_sec0/double(86400);
  // save Matlab time
  t_matlabtime = c_matlabtime;
}


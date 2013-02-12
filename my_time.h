
#ifndef MY_TIME_H
#define MY_TIME_H

// standard headers
#include <string>
#include <ctime>

namespace My {
  // universal class for operations with date and time, much more general and
  // handy than ctime structures
  class Time {
      int _year, _month, _day, _hour, _minute;
      double _second; // seconds including milliseconds
      int _doy; // day of year
      double _unixtime;
      int _local_utc_shift; // shift between local and UTC time in hours
    public:
      // constructors
      Time(); // initialize object with current time
      Time(int, int, int, int, int, double); // initialize object with year, month,
                                          // day, hour, minute and second
      Time(int year, int month, int day, int hour, int minute, int second) {
        Time(year, month, day, hour, minute, double(second));
      }
      Time(int, int, int, int, double); // initialize object with year,
                                     // day of year, hour, minute and second
      Time(int year, int doy, int hour, int minute, int second) {
        Time(year, doy, hour, minute, double(second));
      }
      Time(std::string); // initialize object with time string, the default
                         // format is yyyy-mm-dd HH:MM:SS
      Time(double); // initialize object with Unix timestamp
      // other methods
      double unixtime() const {return _unixtime;} // return Unix time
      int year() const {return _year;} // return current year
      int month() const {return _month;} // return current month
      int day() const {return _day;} // return current day
      int hour() const {return _hour;} // return current hour
      int minute() const {return _minute;} // return current minute
      double second() const {return _second;} // return current second
      int doy() const {return _doy;} // return current day of year
      std::string strtime(); // return formatted time string
      // returns difference between two times in seconds
      double operator-(Time that) const {return _unixtime-that.unixtime();}
      // comparison operators
      bool operator<(Time that) const {return _unixtime < that.unixtime();}
      bool operator>(Time that) const {return _unixtime > that.unixtime();}
      bool operator<=(Time that) const {return _unixtime <= that.unixtime();}
      bool operator>=(Time that) const {return _unixtime >= that.unixtime();}
      Time& add(int, std::string); // add some time to the Time object
    private:
      void initLocalUtcShift(); // initialize shift between local and UTC time
      void initByUnixtime(double); // fill Time object member values using Unix
                                   // time
  };
}
#endif


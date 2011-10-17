
#ifndef MY_TIME_H
#define MY_TIME_H

// standard headers
#include <string>
#include <ctime>

namespace My {
  // universal class for operations with date and time, much more general and
  // handy than ctime structures
  class Time {
      int _year, _month, _day, _hour, _minute, _second;
      int _doy; // day of year
      int _unixtime;
      double _matlabtime;
      int _local_utc_shift; // shift between local and UTC time in hours
    public:
      // constructors
      Time(); // initialize object with current time
      Time(int, int, int, int, int, int); // initialize object with year, month,
                                          // day, hour, minute and second
      Time(int, int, int, int, int); // initialize object with year,
                                     // day of year, hour, minute and second
      Time(std::string); // initialize object with time string, the default
                         // format is yyyy-mm-dd HH:MM:SS
      Time(int, std::string); // initialize object with Unix timestamp
      Time(double, std::string); // initialize object with Matlab timestamp
      // other methods
      int unixtime() const {return _unixtime;} // return Unix time
      double matlabtime() const {return _matlabtime;} // return Matlab time
      int year() const {return _year;} // return current year
      int month() const {return _month;} // return current month
      int day() const {return _day;} // return current day
      int hour() const {return _hour;} // return current hour
      int minute() const {return _minute;} // return current minute
      int second() const {return _second;} // return current second
      int doy() const {return _doy;} // return current day of year
      // returns difference between two times in seconds
      int operator-(Time that) const {return _unixtime-that.unixtime();}
      // comparison operators
      bool operator<(Time that) const {return _unixtime < that.unixtime();}
      bool operator>(Time that) const {return _unixtime > that.unixtime();}
      bool operator<=(Time that) const {return _unixtime <= that.unixtime();}
      bool operator>=(Time that) const {return _unixtime >= that.unixtime();}
      Time& add(int, std::string); // add some time to the Time object
    private:
      void initLocalUtcShift(); // initialize shift between local and UTC time
      void initByUnixtime(time_t); // fill Time object member values using Unix
                                   // time
  };
}
#endif


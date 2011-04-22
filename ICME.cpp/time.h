
#ifndef TIME_H
#define TIME_H

#include <string>
#include <ctime>

using namespace std;

// universal class for operations with date and time, much more general and
// handy than ctyme structures
class Time {
    int t_year, t_month, t_day, t_hour, t_minute, t_second;
    int t_doy; // day of year
    int t_unixtime, t_mysqltime;
    double t_matlabtime;
    string t_stringtime;
    int t_local_utc_shift; // shift between local and UTC time in hours
    tm c_tm;
    time_t c_unixtime;
  public:
    // constructors
    Time(); // initialize object with current time
    Time(int, int, int, int, int, int); // initialize object with year, month,
                                        // day, hour, minute and second
    Time(int, int, int, int, int); // initialize object with year,
                                   // day of year, hour, minute and second

    Time(string); // initialize object with time string, the default format is
                  // yyyy-mm-dd HH:MM:SS

    Time(int, string); // initialize object with Unix, Matlab or
                       // MySQL timestamp
    // other methods
    int unixtime() const {return t_unixtime;} // return unixtime
    double matlabtime() const {return t_matlabtime;} // return Matlab time
    int mysqltime() const {return t_mysqltime;} // return MySQL time
    string stringtime() const {return t_stringtime;} // return formatted
                                                     // string time
    int year() const {return t_year;} // return current year
    int month() const {return t_month;} // return current month
    int day() const {return t_day;} // return current day
    int hour() const {return t_hour;} // return current hour
    int minute() const {return t_minute;} // return current minute
    int second() const {return t_second;} // return current second
    int doy() const {return t_doy;} // return current day of year
    int operator-(Time); // returns difference between two times in seconds
    // comparison operators
    bool operator<(Time);
    bool operator>(Time);
    bool operator<=(Time);
    bool operator>=(Time);
    Time& add(int, string); // add some time to the Time object
  private:
    void initLocalUtcShift(); // initialize shift between local and UTC time
    void initTm(); // initialize ctime strucure
    // calculate Matlab time from ctime strucure and Unix time,
    // service only method
    void calculateMatlabTime();
};

#endif


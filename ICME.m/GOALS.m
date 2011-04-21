%{
Goals:
1. save everything in files, all data, all.
2. use fortran for slow tasks (read data from file, filter it and save to another file for matlab)

Classes and methods

Data::readData(path, timeBegin, timeEnd, gap, )
%}

Config::init()
Config::read()
Config::filter()

Data::init()
Data::read() use C++
Data::filter()
Data::plot()
Data::save()

DhtFrame::init()
DhtFrame::compute() use C++
DhtFrame::plot()
DhtFrame::save()

GsrAxes::init()
GsrAxes::compute() use C++
GsrAxes::plot()
GsrAxes::save()

GsrCurve::init()
GsrCurve::fit()
GsrCurve::plot()
GsrCurve::save()

GsrMap::init()
GsrMap::compute()
GsrMap::plot()
GsrMap::save()

Code all time consuming staff in mex files, use C++ for that.

Eigen+GSL

Where do I actually need a vector container?
1. when reading data from file
2. ?
What maths do I need?
1. all matrix algebra (eigen)
2. eigen values and vectors (eigen)
3. sorting (gsl)
4. polynomial fit
5. minimum and maximum search (eigen, gsl)
6. numerical differentiation
7. numerical integration
8. vector projection (eigen + me)
9. smoothing
10. resampling (gsl)
11. non-linear least squares curve fitting

Directory structure

IcmeLib

IcmeApp
    icme.cpp


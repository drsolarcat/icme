
PROGRAM = icme
BENCHMARK = bench

CXXSOURCES = config.cpp data.cpp my_time.cpp event.cpp \
	mva_analyzer.cpp dht_analyzer.cpp gsr_analyzer.cpp gsr_curve.cpp \
	curve.cpp branched_curve.cpp integrator.cpp differentiator.cpp \
	gsl_fit_epe.cpp gsl_fit_poly.cpp gsl_fit_exp.cpp gnuplot.cpp plotter.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)

CXX = g++
GSL = -lgsl -lgslcblas -lm
BOOST = -lutil -lboost_iostreams
MATLAB = -lmx -leng -lmat -lut
FLAGS = -O3

all: $(PROGRAM) $(BENCHMARK)

$(PROGRAM): $(CXXOBJECTS) $(PROGRAM).o
	$(CXX) -o $@ $@.o $(CXXOBJECTS) $(GSL) $(BOOST) $(MATLAB) $(FLAGS)

$(BENCHMARK): $(CXXOBJECTS) $(BENCHMARK).o
	$(CXX) -o $@ $@.o $(CXXOBJECTS) $(GSL) $(BOOST) $(MATLAB) $(FLAGS)

$(PROGRAM).o: icme.cpp
	$(CXX) -c -o icme.o icme.cpp $(FLAGS)

$(BENCHMARK).o: bench.cpp
	$(CXX) -c -o bench.o bench.cpp $(FLAGS)

config.o: config.h config.cpp
	$(CXX) -c -o config.o config.cpp $(FLAGS)

data.o: data.h data.cpp
	$(CXX) -c -o data.o data.cpp $(FLAGS)

my_time.o: my_time.h my_time.cpp
	$(CXX) -c -o my_time.o my_time.cpp $(FLAGS)

mva_analyzer.o: mva_analyzer.h mva_analyzer.cpp
	$(CXX) -c -o mva_analyzer.o mva_analyzer.cpp $(FLAGS)

dht_analyzer.o: dht_analyzer.h dht_analyzer.cpp
	$(CXX) -c -o dht_analyzer.o dht_analyzer.cpp $(FLAGS)

gsr_analyzer.o: gsr_analyzer.h gsr_analyzer.cpp
	$(CXX) -c -o gsr_analyzer.o gsr_analyzer.cpp $(FLAGS)

curve.o: curve.h curve.cpp
	$(CXX) -c -o curve.o curve.cpp $(FLAGS)

gsr_curve.o: gsr_curve.h gsr_curve.cpp
	$(CXX) -c -o gsr_curve.o gsr_curve.cpp $(FLAGS)

branched_curve.o: branched_curve.h branched_curve.cpp
	$(CXX) -c -o branched_curve.o branched_curve.cpp $(FLAGS)

integrator.o: integrator.h integrator.cpp
	$(CXX) -c -o integrator.o integrator.cpp $(FLAGS)

differentiator.o: differentiator.h differentiator.cpp
	$(CXX) -c -o differentiator.o differentiator.cpp $(FLAGS)

plotter.o: plotter.h plotter.cpp
	$(CXX) -c -o plotter.o plotter.cpp $(FLAGS)

gsl_fit_epe.o: gsl_fit_epe.h gsl_fit_epe.cpp
	$(CXX) -c -o gsl_fit_epe.o gsl_fit_epe.cpp $(FLAGS)

gsl_fit_poly.o: gsl_fit_poly.h gsl_fit_poly.cpp
	$(CXX) -c -o gsl_fit_poly.o gsl_fit_poly.cpp $(FLAGS)

gsl_fit_exp.o: gsl_fit_exp.h gsl_fit_exp.cpp
	$(CXX) -c -o gsl_fit_exp.o gsl_fit_exp.cpp $(FLAGS)

gnuplot.o: gnuplot.h gnuplot.cpp
	$(CXX) -c -o gnuplot.o gnuplot.cpp $(FLAGS)

event.o: event.h event.cpp
	$(CXX) -c -o event.o event.cpp $(FLAGS)

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM).o $(BENCHMARK).o $(PROGRAM) $(BENCHMARK)

run:
	./$(PROGRAM)

test:
	./$(BENCHMARK)



PROGRAM = icme
BENCHMARK = bench

CXXSOURCES = config.cpp data.cpp my_time.cpp event.cpp \
	mva_analyzer.cpp dht_analyzer.cpp gsr_analyzer.cpp gsr_curve.cpp \
	curve.cpp branched_curve.cpp integrator.cpp differentiator.cpp \
	fit_poly.cpp fit_exp.cpp fit_cexp.cpp fit_polyexp.cpp gnuplot.cpp \
	plotter.cpp
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

fit_poly.o: fit_poly.h fit_poly.cpp
	$(CXX) -c -o fit_poly.o fit_poly.cpp $(FLAGS)

fit_exp.o: fit_exp.h fit_exp.cpp
	$(CXX) -c -o fit_exp.o fit_exp.cpp $(FLAGS)

fit_cexp.o: fit_cexp.h fit_cexp.cpp
	$(CXX) -c -o fit_cexp.o fit_cexp.cpp $(FLAGS)

fit_polyexp.o: fit_polyexp.h fit_polyexp.cpp
	$(CXX) -c -o fit_polyexp.o fit_polyexp.cpp $(FLAGS)

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



# check where we are compiling the project
ifeq ($(USERNAME),isavnin) # at work
	CPLUS_INCLUDE_PATH=/home/isavnin/usr/local/include:/usr/local/matlab2009b/extern/include
	LIBRARY_PATH=/home/isavnin/usr/local/lib:/usr/local/matlab2009b/extern/lib/glnxa64:/usr/local/matlab2009b/bin/glnxa64
	LD_LIBRARY_PATH=/home/isavnin/usr/local/lib:/usr/local/matlab2009b/bin/glnxa64
endif

#ifeq ($(USERNAME),inferno) # on my desktop
#	CPLUS_INCLUDE_PATH=/home/inferno/usr/local/include:/usr/local/MATLAB/R2010b/extern/include
#	LIBRARY_PATH=/home/inferno/usr/local/lib:/usr/local/MATLAB/R2010b/extern/lib/glnxa64:/usr/local/MATLAB/R2010b/bin/glnxa64
#	LD_LIBRARY_PATH=/home/inferno/usr/local/lib:/usr/local/MATLAB/R2010b/bin/glnxa64
#endif

ifeq ($(USERNAME),inferno) # on my laptop
	CPLUS_INCLUDE_PATH=/home/inferno/usr/local/include:/usr/local/MATLAB/R2011a/extern/include
	LIBRARY_PATH=/home/inferno/usr/local/lib:/usr/local/MATLAB/R2011a/extern/lib/glnx86:/usr/local/MATLAB/R2011a/bin/glnx86
	LD_LIBRARY_PATH=/home/inferno/usr/local/lib:/usr/local/MATLAB/R2011a/bin/glnx86
endif

# export the paths
export CPLUS_INCLUDE_PATH LIBRARY_PATH

# names of executables
PROGRAM = icme
BENCHMARK = bench

# sources
CXXSOURCES = config.cpp data.cpp my_time.cpp event.cpp \
	mva_analyzer.cpp dht_analyzer.cpp gsr_analyzer.cpp gsr_curve.cpp \
	curve.cpp branched_curve.cpp integrator.cpp differentiator.cpp \
	fit_abstract.cpp fit_poly.cpp fit_exp.cpp fit_poly_exp.cpp gnuplot.cpp \
	plotter.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)

#flags
CXX = g++
GSL = -lgsl -lgslcblas -lm
BOOST = -lutil -lboost_iostreams
MATLAB = -lmx -leng -lmat -lut
LOG4CPLUS = -llog4cplus
CXXFLAGS = -O3 -Wl,-rpath,$(LD_LIBRARY_PATH),-rpath-link,$(LD_LIBRARY_PATH)

all: $(PROGRAM) $(BENCHMARK)

$(PROGRAM): $(CXXOBJECTS) $(PROGRAM).o
	$(CXX) -o $@ $@.o $(CXXOBJECTS) $(GSL) $(BOOST) $(MATLAB) $(LOG4CPLUS) $(CXXFLAGS)

$(BENCHMARK): $(CXXOBJECTS) $(BENCHMARK).o
	$(CXX) -o $@ $@.o $(CXXOBJECTS) $(GSL) $(BOOST) $(MATLAB) $(LOG4CPLUS) $(CXXFLAGS)

$(PROGRAM).o: icme.cpp
	$(CXX) -c -o icme.o icme.cpp $(CXXFLAGS)

$(BENCHMARK).o: bench.cpp
	$(CXX) -c -o bench.o bench.cpp $(CXXFLAGS)

config.o: config.h config.cpp
	$(CXX) -c -o config.o config.cpp $(CXXFLAGS)

data.o: data.h data.cpp
	$(CXX) -c -o data.o data.cpp $(CXXFLAGS)

my_time.o: my_time.h my_time.cpp
	$(CXX) -c -o my_time.o my_time.cpp $(CXXFLAGS)

mva_analyzer.o: mva_analyzer.h mva_analyzer.cpp
	$(CXX) -c -o mva_analyzer.o mva_analyzer.cpp $(CXXFLAGS)

dht_analyzer.o: dht_analyzer.h dht_analyzer.cpp
	$(CXX) -c -o dht_analyzer.o dht_analyzer.cpp $(CXXFLAGS)

gsr_analyzer.o: gsr_analyzer.h gsr_analyzer.cpp
	$(CXX) -c -o gsr_analyzer.o gsr_analyzer.cpp $(CXXFLAGS)

curve.o: curve.h curve.cpp
	$(CXX) -c -o curve.o curve.cpp $(CXXFLAGS)

gsr_curve.o: gsr_curve.h gsr_curve.cpp
	$(CXX) -c -o gsr_curve.o gsr_curve.cpp $(CXXFLAGS)

branched_curve.o: branched_curve.h branched_curve.cpp
	$(CXX) -c -o branched_curve.o branched_curve.cpp $(CXXFLAGS)

integrator.o: integrator.h integrator.cpp
	$(CXX) -c -o integrator.o integrator.cpp $(CXXFLAGS)

differentiator.o: differentiator.h differentiator.cpp
	$(CXX) -c -o differentiator.o differentiator.cpp $(CXXFLAGS)

plotter.o: plotter.h plotter.cpp
	$(CXX) -c -o plotter.o plotter.cpp $(CXXFLAGS)

fit_abstract.o: fit_abstract.h fit_abstract.cpp
	$(CXX) -c -o fit_abstract.o fit_abstract.cpp $(CXXFLAGS)

fit_poly.o: fit_poly.h fit_poly.cpp
	$(CXX) -c -o fit_poly.o fit_poly.cpp $(CXXFLAGS)

fit_exp.o: fit_exp.h fit_exp.cpp
	$(CXX) -c -o fit_exp.o fit_exp.cpp $(CXXFLAGS)

fit_poly_exp.o: fit_poly_exp.h fit_poly_exp.cpp
	$(CXX) -c -o fit_poly_exp.o fit_poly_exp.cpp $(CXXFLAGS)

gnuplot.o: gnuplot.h gnuplot.cpp
	$(CXX) -c -o gnuplot.o gnuplot.cpp $(CXXFLAGS)

event.o: event.h event.cpp
	$(CXX) -c -o event.o event.cpp $(CXXFLAGS)

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM).o $(BENCHMARK).o $(PROGRAM) $(BENCHMARK)

run:
	./$(PROGRAM)

test:
	./$(BENCHMARK)


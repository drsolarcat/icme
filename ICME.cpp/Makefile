
PROGRAM = icme

CXXSOURCES = icme.cpp config.cpp data.cpp time.cpp event.cpp \
	mva_analyzer.cpp dht_analyzer.cpp gsr_analyzer.cpp gsr_curve.cpp \
	curve.cpp branched_curve.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)
CXX = g++
GSL = -lgsl -lgslcblas -lm

all: $(PROGRAM)

$(PROGRAM): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS) $(GSL)

config.o: config.h config.cpp
	$(CXX) -c -o config.o config.cpp

data.o: data.h data.cpp
	$(CXX) -c -o data.o data.cpp

time.o: time.h time.cpp
	$(CXX) -c -o time.o time.cpp

event.o: event.h event.cpp
	$(CXX) -c -o event.o event.cpp

mva_analyzer.o: mva_analyzer.h mva_analyzer.cpp
	$(CXX) -c -o mva_analyzer.o mva_analyzer.cpp

dht_analyzer.o: dht_analyzer.h dht_analyzer.cpp
	$(CXX) -c -o dht_analyzer.o dht_analyzer.cpp

gsr_analyzer.o: gsr_analyzer.h gsr_analyzer.cpp
	$(CXX) -c -o gsr_analyzer.o gsr_analyzer.cpp

curve.o: curve.h curve.cpp
	$(CXX) -c -o curve.o curve.cpp

gsr_curve.o: gsr_curve.h gsr_curve.cpp
	$(CXX) -c -o gsr_curve.o gsr_curve.cpp

branched_curve.o: branched_curve.h branched_curve.cpp
	$(CXX) -c -o branched_curve.o branched_curve.cpp

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM)

run:
	./$(PROGRAM)


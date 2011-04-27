
PROGRAM = icme

CXXSOURCES = icme.cpp config.cpp data.cpp time.cpp event.cpp mva_analyzer.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)
CXX = g++

all: $(PROGRAM)

$(PROGRAM): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS)

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

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM)

run:
	./$(PROGRAM)


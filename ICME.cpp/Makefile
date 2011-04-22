
PROGRAM = icme

CXXSOURCES = icme.cpp config.cpp data.cpp time.cpp
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

time_test: time_test.cpp time.o
	$(CXX) -o time_test time_test.cpp time.o

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM)

run:
	./$(PROGRAM)


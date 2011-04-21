
PROGRAM = icme

CXXSOURCES = icme.cpp config.cpp
CXXOBJECTS = $(CXXSOURCES:.cpp=.o)
CXX = g++

all: $(PROGRAM)

$(PROGRAM): $(CXXOBJECTS)
	$(CXX) -o $@ $(CXXOBJECTS)

config.o: config.h config.cpp
	$(CXX) -c -o config.o config.cpp

clean:
	$(RM) -f $(CXXOBJECTS) $(PROGRAM)

run:
	./$(PROGRAM)


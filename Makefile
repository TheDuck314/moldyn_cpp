OBJS = $(addsuffix .o, $(basename $(wildcard *.cpp)))

CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -fopenmp
#CXXFLAGS = -pg -std=c++11 -Wall -fopenmp
BINARY = moldyn

all: $(OBJS)
	g++ $(CXXFLAGS) -o $(BINARY) $(OBJS)

clean:
	rm -f $(OBJS)
	rm $(BINARY)

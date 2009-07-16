CXX = g++
CXX_FLAGS = -O3 -arch x86_64 -Wall -lpthread

all : pointpacker gen_boundary clean

pointpacker : geometry.o pointpacker.o
	$(CXX) pointpacker.o geometry.o -o pointpacker $(CXX_FLAGS)

gen_boundary: geometry.o gen_boundary.o
	$(CXX) gen_boundary.o geometry.o -o gen_boundary $(CXX_FLAGS)

pointpacker.o : pointpacker.cpp geometry.h
	$(CXX) -c pointpacker.cpp $(CXX_FLAGS)

geometry.o : geometry.cpp geometry.h
	$(CXX) -c geometry.cpp $(CXX_FLAGS)

gen_boundary.o : geometry.h gen_boundary.cpp
	$(CXX) -c gen_boundary.cpp $(CXX_FLAGS)

.PHONY: clean
clean :
	-rm -f *.o *~; true
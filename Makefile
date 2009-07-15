CXX = g++
CXX_FLAGS = -O3 -arch x86_64 -Wall -lpthread

.PHONY: all
all : pointpacker clean

pointpacker : geometry.o pointpacker.o
	$(CXX) pointpacker.o geometry.o -o pointpacker $(CXX_FLAGS)

pointpacker.o : pointpacker.cpp geometry.h
	$(CXX) -c pointpacker.cpp $(CXX_FLAGS)

geometry.o : geometry.cpp geometry.h
	$(CXX) -c geometry.cpp $(CXX_FLAGS)

.PHONY: clean
clean :
	-rm -f *.o *~; true
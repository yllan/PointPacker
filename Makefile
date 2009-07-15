CXX = g++
CXX_FLAGS = -O3 -arch x86_64 -Wall

.PHONY: all
all : bruteforce clean

bruteforce : geometry.o bruteforce.o
	$(CXX) bruteforce.o geometry.o -o bruteforce $(CXX_FLAGS)

bruteforce.o : bruteforce.cpp geometry.h
	$(CXX) -c bruteforce.cpp $(CXX_FLAGS)

geometry.o : geometry.cpp geometry.h
	$(CXX) -c geometry.cpp $(CXX_FLAGS)

.PHONY: clean
clean :
	-rm -f *.o *~; true
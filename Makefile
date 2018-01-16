CXXFLAGS=-lwat -I$(HOME)/soft/eigen-3.3.4/ -O2 -pg
test:posimp.o test.o
	g++ posimp.o test.o -o test $(CXXFLAGS)
posimp.o:posimp.cpp posimp.hpp
test.o:test.cpp posimp.hpp

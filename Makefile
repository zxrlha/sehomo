CXXFLAGS=-lwat -I$(HOME)/soft/eigen-3.3.4/ -Ofast -mtune=native -march=native -g -I$(HOME)/projects/atohoc/athl/src -g -I/usr/include/flint -I/usr/include/arb -Wall -pg
LDFLAGS=-L$(HOME)/projects/atohoc/athl/src/.libs -lathl -lstdc++ -lm -lflint -larb -lgmp -lmpfr -pg
test_BG:test_BG.o BG.o
test:posimp.o test.o lastimp.o phyimp.o
	g++ posimp.o test.o lastimp.o phyimp.o -o test $(CXXFLAGS) -lboost_serialization
posimp.o:posimp.cpp posimp.hpp
lastimp.o:lastimp.cpp posimp.hpp lastimp.hpp
test.o:test.cpp posimp.hpp lastimp.hpp finalstate.hpp phyimp.hpp
phyimp.o:phyimp.hpp phyimp.cpp

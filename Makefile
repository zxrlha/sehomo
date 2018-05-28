CXXFLAGS=`pkg-config eigen3 --cflags` -O2 -mtune=native -march=native -g
all:prep_ic run_phy
prep_ic:posimp.o prep_ic.o comimp.o timer.o
	g++ -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -lboost_serialization
run_phy:phyimp.o run_phy.o timer.o
	g++ -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -lboost_serialization
posimp.o:posimp.cpp posimp.hpp timer.hpp
phyimp.o:phyimp.cpp posimp.hpp phyimp.hpp timer.hpp
comimp.o:comimp.hpp comimp.cpp
prep_ic.o:prep_ic.cpp rcm.hpp lorentz.hpp

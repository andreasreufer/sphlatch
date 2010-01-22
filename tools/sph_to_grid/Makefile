all: sphtogrid

sphtogrid:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
	  -I../../src/ -fopenmp \
	  -o sph2grid sph2grid.cpp


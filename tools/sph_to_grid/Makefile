all: sphtogrid_CSS sphtogrid_CSV sphtogrid_MSS sphtogrid_MSV

sphtogrid_CSS:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
	  -I../../src/ -fopenmp \
	  -DSPHLATCH_INTERPOLATE_SCALAR \
	  -o sph2grid_CSS sph2grid.cpp

sphtogrid_CSV:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
	  -I../../src/ -fopenmp \
	  -DSPHLATCH_INTERPOLATE_VECTOR \
	  -o sph2grid_CSV sph2grid.cpp

sphtogrid_MSS:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
	  -I../../src/ -fopenmp \
	  -DSPHLATCH_MISCIBLE \
	  -DSPHLATCH_INTERPOLATE_SCALAR \
	  -o sph2grid_MSS sph2grid.cpp

sphtogrid_MSV:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
	  -I../../src/ -fopenmp \
	  -DSPHLATCH_MISCIBLE \
	  -DSPHLATCH_INTERPOLATE_VECTOR \
	  -o sph2grid_MSV sph2grid.cpp


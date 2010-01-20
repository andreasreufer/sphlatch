all:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src \
	  -fopenmp -DSPHLATCH_OPENMP\
	  -lhdf5 -lz \
	  -o hdf5simple hdf5simple.cpp

clean:
	rm hdf5simple 2>&1 >/dev/null; true

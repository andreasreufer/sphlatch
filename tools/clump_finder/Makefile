all: find_clumps

find_clumps:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -o find_clumps find_clumps.cpp

clean:
	rm find_clumps 2>&1 >/dev/null; true


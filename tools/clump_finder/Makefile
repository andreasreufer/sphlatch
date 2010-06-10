all: find_clumps_friends find_clumps_potential

find_clumps_friends:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DCLUMPFIND_FOF \
	   -o find_clumps_friends find_clumps.cpp

find_clumps_potential:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -o find_clumps_potential find_clumps.cpp

clean:
	rm find_clumps_friends find_clumps_potential 2>&1 >/dev/null; true


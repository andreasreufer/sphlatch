all: clump_finder_friends clump_finder_potential

clump_finder_friends:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DCLUMPFIND_FOF \
	   -o clump_finder_friends clump_finder.cpp

clump_finder_potential:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -o clump_finder_potential clump_finder.cpp

clean:
	rm clump_finder_friends clump_finder_potential 2>&1 >/dev/null; true


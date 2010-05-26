all: clump_finder_friends

clump_finder_friends:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -o clump_finder_friends clump_finder_friends.cpp

clean:
	rm clump_finder_friends 2>&1 >/dev/null; true


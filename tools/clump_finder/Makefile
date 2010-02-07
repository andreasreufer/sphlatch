all: clump_finder

clump_finder:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -o clump_finder clump_finder.cpp

clean:
	rm clump_finder 2>&1 >/dev/null; true


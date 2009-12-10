all: bfnbody

bfnbody:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -o bf_nbody bf_nbody.cpp
	
	
clean:
	rm bf_nbody 2>&1 >/dev/null; true


all: treefun

treefun:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -lhdf5 -lz \
	  -fopenmp \
	  -DSPHLATCH_GRAVITY \
	  -o treeFun treeFun.cpp

nodesfun:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -o nodesFun nodesFun.cpp

clean:
	rm treeFun nodesFun 2>&1 >/dev/null; true

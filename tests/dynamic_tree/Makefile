all: treefun

treefun:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -lhdf5 -lz -fopenmp\
	  -DSPHLATCH_GRAVITY \
	  -o treeFun treeFun.cpp

clean:
	rm treeFun 2>&1 >/dev/null; true

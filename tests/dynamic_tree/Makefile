all: treefun nodesize buildcz

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

nodesize:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -fopenmp \
	  -o nodeSize nodeSize.cpp

buildcz:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -fopenmp \
	  -lhdf5 -lz \
	  -o buildCZ buildCZ.cpp

workertest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -fopenmp \
	  -lhdf5 -lz \
	  -o workerTest workerTest.cpp

clean:
	rm treeFun nodesFun 2>&1 >/dev/null; true

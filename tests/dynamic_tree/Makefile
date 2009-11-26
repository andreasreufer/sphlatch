all: treefun nodesize buildcz nodesfun workertest

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


bhtree:
	        $(CXX) $(CXXFLAGS) $(LDFLAGS) \
		  -I../../src -I../../src/dynamic_tree \
		  -lhdf5 -lz \
		  -fopenmp \
		  -o treeFun treeFun.cpp

sphworker:
	        $(CXX) $(CXXFLAGS) $(LDFLAGS) \
		  -I../../src -I../../src/dynamic_tree \
		  -lhdf5 -lz \
		  -fopenmp \
		  -o sphWorker sphWorker.cpp

particleset:
	        $(CXX) $(CXXFLAGS) $(LDFLAGS) \
		  -I../../src -I../../src/dynamic_tree \
		  -lhdf5 -lz \
		  -fopenmp \
		  -o particleSet particleSet.cpp


includefun:
	        $(CXX) $(CXXFLAGS) $(LDFLAGS) \
		  -I../../src -I../../src/dynamic_tree \
		  -lhdf5 -lz \
		  -fopenmp \
		  -o includeFun include.cpp
		
nBody:
	        $(CXX) $(CXXFLAGS) $(LDFLAGS) \
		  -I../../src -I../../src/dynamic_tree \
		  -lhdf5 -lz \
		  -fopenmp \
		  -o NBody nbody.cpp


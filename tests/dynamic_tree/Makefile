all: nodesfun nodesize workertest sphworker particleset nbody

nodesfun:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src \
	  -o nodesFun nodesFun.cpp

nodesize:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src \
	  -fopenmp \
	  -o nodeSize nodeSize.cpp

workertest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src \
	  -fopenmp \
	  -lhdf5 -lz \
	  -o workerTest workerTest.cpp

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


nbody:
	        $(CXX) $(CXXFLAGS) $(LDFLAGS) \
		  -I../../src \
		  -lhdf5 -lz \
		  -fopenmp \
		  -o nBody nbody.cpp

clean:
	rm nodesFun nodeSize workerTest sphWorker \
	  particleSet nBody 2>&1 >/dev/null; true


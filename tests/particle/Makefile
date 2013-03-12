all: parttest

parttest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -I../../src -I../../src/dynamic_tree \
	  -fopenmp \
	  -lhdf5 -lz \
	  -o particleTest particleTest.cpp

clean:
	rm particleTest 2>&1 >/dev/null; true

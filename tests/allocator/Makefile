all: allocatorTest

allocatorTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
	  -I../../src/dynamic_tree \
	  -lhdf5 -lz\
	  -DSPHLATCH_GRAVITY \
	  -o allocator_test allocator_test.cpp


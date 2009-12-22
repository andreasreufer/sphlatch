all: allocatorTest

allocatorTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
	  -lhdf5 -lz\
	  -DSPHLATCH_GRAVITY \
	  -o allocator_test allocator_test.cpp


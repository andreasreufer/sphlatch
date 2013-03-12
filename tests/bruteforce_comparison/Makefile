all: bfcompS bfcomp_

bfcompS:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -DSPHLATCH_GRAVITY \
	  -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	  -fopenmp \
	  -o bfcomp_S bfcomp.cpp

bfcomp_:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -DSPHLATCH_GRAVITY \
	  -fopenmp \
	  -o bfcomp__ bfcomp.cpp


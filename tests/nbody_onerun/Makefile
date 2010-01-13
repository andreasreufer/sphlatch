all: onerun_S onerun__

onerun_S:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -DSPHLATCH_GRAVITY \
	  -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	  -fopenmp \
	  -o nbody_onerun_S nbody_onerun.cpp
onerun__:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
	  -lhdf5 -lz -lboost_program_options \
	  -DSPHLATCH_GRAVITY \
	  -DSPHLATCH_GRAVITY_NOSMOOTHING \
	  -fopenmp \
	  -o nbody_onerun__ nbody_onerun.cpp


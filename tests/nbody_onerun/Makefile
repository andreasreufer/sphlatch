all: onerun_S onerun__

onerun_S:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
	  -lhdf5 -lz -lboost_program_options \
	  -DSPHLATCH_GRAVITY \
	  -DSPHLATCH_GRAVSPLINESMOOTH \
	  -o nbody_onerun_S nbody_onerun.cpp
onerun__:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
	  -lhdf5 -lz -lboost_program_options \
	  -DSPHLATCH_GRAVITY \
	  -o nbody_onerun__ nbody_onerun.cpp


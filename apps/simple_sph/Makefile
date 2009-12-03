all: simple_sph_______ simple_sph_GHU___

simple_sph_______:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -o simple_sph_______ simple_sph.cpp

simple_sph_GHU___:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -fopenmp \
	   -o simple_sph_______ simple_sph.cpp


clean:
	rm simple_sph_* 2>&1 >/dev/null; true


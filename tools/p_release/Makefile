all: getPressures

getPressures:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_ANEOS \
	   -o p_release p_release.cpp
clean:
	rm p_release 2>&1 >/dev/null; true


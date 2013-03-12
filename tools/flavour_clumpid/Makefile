all: flavourClumpID

flavourClumpID:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_ANEOS \
	   -o flavour_clumpid flavour_clumpid.cpp
clean:
	rm p_release 2>&1 >/dev/null; true


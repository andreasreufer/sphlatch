all: moonPostProc

moonPostProc:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src -fopenmp \
				      -DSPHLATCH_ANEOS \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_POTENTIAL\
				      -o moon_pp moon_postproc.cpp

clean:
	rm moon_pp 2>&1 >/dev/null;\
	true

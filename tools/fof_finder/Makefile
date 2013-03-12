all: findFOF_A_

findFOF_A_:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_ANEOS_TRACKORIG \
	   -o find_fof find_fof.cpp

clean:
	rm find_clumps_?? 2>&1 >/dev/null; true


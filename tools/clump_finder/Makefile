all: findClumps_A_

findClumps___:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -o find_clumps__ find_clumps.cpp

findClumps_AG:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_CALC_POTENTIAL \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_ANEOS \
	   -o find_clumps_AG find_clumps.cpp

findClumps_A_:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_FIND_CLUMPS \
	   -DSPHLATCH_ANEOS \
	   -o find_clumps_A_ find_clumps.cpp

clean:
	rm find_clumps_?? 2>&1 >/dev/null; true


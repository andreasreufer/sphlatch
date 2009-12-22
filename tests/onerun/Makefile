all: onerun_GHUA__


onerun_GHUA__: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -fopenmp \
	   -o onerun_GHUA__ onerun.cpp
libaneos:
	    $(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o


clean:
	rm onerun_* 2>&1 >/dev/null; true


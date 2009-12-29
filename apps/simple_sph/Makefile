all: simple_sph_______ simple_sph_GHU___ simple_sph_GHUA__ simple_sph_GHUA_F simple_sph_GHUa__ simple_sph_GHUa_F

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
	   -o simple_sph_GHU___ simple_sph.cpp

simple_sph_GHUA__: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -fopenmp \
	   -o simple_sph_GHUA__ simple_sph.cpp

simple_sph_GHUa__: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -fopenmp \
	   -o simple_sph_GHUa__ simple_sph.cpp

simple_sph_GHUA_F: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -fopenmp \
	   -o simple_sph_GHUA_F simple_sph.cpp

simple_sph_GHUa_F: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_FRICTION \
	   -fopenmp \
	   -o simple_sph_GHUa_F simple_sph.cpp


libaneos:
	    $(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o


clean:
	rm simple_sph_* 2>&1 >/dev/null; true


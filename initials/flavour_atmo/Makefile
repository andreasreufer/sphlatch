all: libaneos flavAtmoH flavAtmoM

flavAtmoH:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -fopenmp \
	  -DSPHLATCH_ANEOS_TABLE \
	  -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	  -o flavour_atmo__ flavour_atmo.cpp

flavAtmoM:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -fopenmp \
	  -DSPHLATCH_ANEOS_TABLE \
	  -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	  -DSPHLATCH_MISCIBLE \
	  -o flavour_atmo_M flavour_atmo.cpp

libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

	

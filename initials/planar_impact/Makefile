all: libmaneos
	  $(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
	    -DSPHLATCH_ANEOS \
	    -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
      	    -o planar_impact planar_impact.cpp
libmaneos:
	  $(FC) $(FFLAGS) -c ../../aux/libaneos/libmaneos.f -o libmaneos.o

	

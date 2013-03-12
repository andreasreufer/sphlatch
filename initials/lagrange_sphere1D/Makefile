all: lagrange1D_A
	
libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

lagrange1D_A: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_ANEOS_TABLE \
				      -o lagrange_sphere1D_A \
				      lagrange_sphere1D.cpp

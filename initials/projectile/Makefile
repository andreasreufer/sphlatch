all: tillotson aneos

tillotson:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
				      -lhdf5 -lz \
				      -lboost_program_options \
				      -o projectile projectile.cpp
	
aneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -lhdf5 -lz \
				      -lboost_program_options \
				      -o projectile_aneos projectile.cpp


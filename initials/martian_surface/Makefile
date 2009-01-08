all: tillotson aneos

tillotson:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o martian_surface martian_surface.cpp
	
aneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -o martian_surface_aneos martian_surface.cpp

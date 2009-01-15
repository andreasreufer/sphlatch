all: h5partAddANEOS h5partAddANEOStable

h5partAddANEOS: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o h5part_add_ANEOSout h5part_add_ANEOSout.cpp

h5partAddANEOStable: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_ANEOS_TABLE \
				      -o h5part_add_ANEOSout_table h5part_add_ANEOSout.cpp

libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

clean:
	rm h5part_add_ANEOSout 2>&1 >/dev/null;\
	rm h5part_add_ANEOSout_table 2>&1 >/dev/null;\
	rm libaneos.o 2>&1 >/dev/null | true


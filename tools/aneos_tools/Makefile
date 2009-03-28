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

h5partGetPhaseHistogram:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -o h5part_get_phase_histogram h5part_getPhaseHistogram.cpp

phaseDiagram: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o aneos_phase_diagram phase_diagram.cpp

libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

clean:
	rm h5part_add_ANEOSout 2>&1 >/dev/null;\
	rm h5part_add_ANEOSout_table 2>&1 >/dev/null;\
	rm libaneos.o 2>&1 >/dev/null | true


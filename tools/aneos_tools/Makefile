all: rootpM rootpA rootuM rootuA generateTablesA generateTablesM phaseDiagram


rootpM: libmaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src \
				      -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
				      -o rootPM rootP.cpp
rootpA: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o rootPA rootP.cpp
rootuM: libmaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src \
				      -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
				      -o rootUM rootU.cpp
rootuA: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o rootUA rootU.cpp

rootsM: libmaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src \
				      -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
				      -o rootSM rootS.cpp
rootsA: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz \
				      -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o rootSA rootS.cpp

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

phaseDiagramM: libmaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libmaneos.o \
				      -o Maneos_phase_diagram phase_diagram.cpp

resetTmaxuav:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -I../../src \
				      -o reset_Tmax_uav reset_Tmax_uav.cpp

generateTablesA: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -I../../src \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o generateTablesA generateTables.cpp

generateTablesM: libmaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -I../../src \
				      -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
				      -o generateTablesM generateTables.cpp


libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

libmaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libmaneos.f -o libmaneos.o

clean:
	rm h5part_add_ANEOSout 2>&1 >/dev/null;\
	rm h5part_add_ANEOSout_table 2>&1 >/dev/null;\
	rm libaneos.o 2>&1 >/dev/null | true


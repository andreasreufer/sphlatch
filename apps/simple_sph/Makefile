all: idealgas tillotson shocktube

idealgas: simple_sph_______ simple_sph_G_____ simple_sph_GH____ \
	  simple_sph_GHU___ simple_sph__HU___

tillotson: simple_sph_GHUT__ simple_sph__HUT__ \
	   simple_sph_GHUTR_ simple_sph__HUTR_ \
	   simple_sph_GHUT_F simple_sph_GHUTRF

aneos:     simple_sph_GHUA__ simple_sph__HUA__ simple_sph_GHUA_F


shocktube:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_SHOCKTUBE \
				      -o shocktube simple_sph.cpp
simple_sph_______:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -o simple_sph_______ simple_sph.cpp

simple_sph_G_____:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -o simple_sph_G_____ simple_sph.cpp

simple_sph_GH____:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -o simple_sph_GH____ simple_sph.cpp

simple_sph_GHU___:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -o simple_sph_GHU___ simple_sph.cpp

simple_sph__HU___:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -o simple_sph__HU___ simple_sph.cpp

simple_sph_GHUT__:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -o simple_sph_GHUT__ simple_sph.cpp

simple_sph__HUT__:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -o simple_sph__HUT__ simple_sph.cpp

simple_sph_GHUTR_:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_INTEGRATERHO \
				      -o simple_sph_GHUTR_ simple_sph.cpp

simple_sph__HUTR_:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_INTEGRATERHO \
				      -o simple_sph__HUTR_ simple_sph.cpp

simple_sph_GHUT_F:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_FRICTION \
				      -o simple_sph_GHUT_F simple_sph.cpp

simple_sph_GHUTRF:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_INTEGRATERHO \
				      -DSPHLATCH_FRICTION \
				      -o simple_sph_GHUTRF simple_sph.cpp

simple_sph_GHUA__: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_ANEOS_TABLE \
				      -o simple_sph_GHUA__ simple_sph.cpp

simple_sph__HUA__: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_ANEOS_TABLE \
				      -o simple_sph__HUA__ simple_sph.cpp

simple_sph_GHUA_F: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_ANEOS_TABLE \
				      -DSPHLATCH_FRICTION \
				      -DSPHLATCH_NONEGPRESS \
				      -o simple_sph_GHUA_F simple_sph.cpp

libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

clean:
	rm simple_sph_* 2>&1 >/dev/null;\
	rm libaneos.o 2>&1 >/dev/null;\
	rm shocktube 2>&1 >/dev/null | true


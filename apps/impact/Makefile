all: impact_GHUA__ impact_GHUT__ impact_GHUa_K

paramsParse:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src \
				      -o paramsParse params.cpp

impact_GHUA__: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_REMOVEESCAPING \
				      -DSPHLATCH_NONEGPRESS \
				      -DSPHLATCH_GRAVSPLINESMOOTH \
				      -o impact_GHUA__ impact.cpp
impact_GHUT__:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_REMOVEESCAPING \
				      -DSPHLATCH_GRAVSPLINESMOOTH \
				      -DSPHLATCH_NONEGPRESS \
				      -o impact_GHUT__ impact.cpp

impact_GHUa_K: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -DSPHLATCH_GRAVSPLINESMOOTH \
				      -DSPHLATCH_NONEGPRESS \
				      -o impact_GHUa_K impact.cpp

libaneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

clean:
	rm impact_* 2>&1 >/dev/null;\
	rm paramParse 2>&1 >/dev/null;\
	rm libaneos.o 2>&1 >/dev/null;\
	true

all: simpleSph shockTube

simpleSph:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o simple_sph simple_sph.cpp

shockTube:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o shocktube shocktube.cpp


variants:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -o simple_sph____ simple_sph.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -o simple_sph___G simple_sph.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -o simple_sph__HG simple_sph.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -o simple_sph_UHG simple_sph.cpp


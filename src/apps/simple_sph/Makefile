all: simpleSph shockTube
variants: simpleSph simpleSphG simpleSphGH simpleSphGHU simpleSph_HU

shockTube:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o shocktube shocktube.cpp


simpleSph:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -o simple_sph____ simple_sph.cpp
simpleSphG:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -o simple_sph___G simple_sph.cpp
simpleSphGH:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -o simple_sph__HG simple_sph.cpp
simpleSphGHU:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -o simple_sph_UHG simple_sph.cpp
simpleSph_HU:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -o simple_sph_UH_ simple_sph.cpp

all: variants shockTube tillotson
variants: simpleSph simpleSphG simpleSphGH simpleSphGHU simpleSph_HU
tillotson: tillotson_sumrho tillotson_intrho

shockTube:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_SHOCKTUBE \
				      -o shocktube simple_sph.cpp
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
tillotson_sumrho:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -o simple_sph_tillotson_sumrho simple_sph.cpp

tillotson_intrho:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TIMEDEP_ENERGY \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_INTEGRATERHO \
				      -o simple_sph_tillotson_intrho simple_sph.cpp

tillotson_intrho_isoT:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_TIMEDEP_SMOOTHING \
				      -DSPHLATCH_TILLOTSON \
				      -DSPHLATCH_INTEGRATERHO \
				      -o simple_sph_tillotson_intrho_isoT simple_sph.cpp

clean:
	rm shocktube simple_sph simple_sph____ simple_sph___G simple_sph__HG \
	simple_sph_UHG simple_sph_UH_ simple_sph_tillotson_sumrho \
	simple_sph_tillotson_intrho simple_sph_tillotson_intrho_isoT \
	2>&1 >/dev/null | true

all: nbody

nbody:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 \
				      -lz -lboost_program_options \
				      -DSPHLATCH_GRAVITY \
				      -DSPHLATCH_GRAVITY_EPSSMOOTHING \
				      -o nbody_integrator nbody_integrator.cpp
	



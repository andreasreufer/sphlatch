all: getEnergies___ getEnergies__G

getEnergies___:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -o get_energies__ get_energies.cpp

getEnergies__G:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_CALC_POTENTIAL \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_ANEOS \
	   -o get_energies_G get_energies.cpp
clean:
	rm get_energies_* 2>&1 >/dev/null; true


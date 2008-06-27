all: simplesph ShockTube

simplesph:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o simple_sph simple_sph.cpp

ShockTube:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o shocktube shocktube.cpp



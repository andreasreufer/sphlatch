all: flavH flavM

flavH:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -fopenmp \
	  -o flavour_sphere_H flavour_sphere.cpp
	
flavM:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -lz -I../../src \
	  -fopenmp \
	  -DSPHLATCH_MISCIBLE \
	  -o flavour_sphere_M flavour_sphere.cpp
	

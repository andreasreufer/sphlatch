all: noDensity

noDensity:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -o no_density no_density.cpp

tomasciiToh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -o tomascii2h5part tomascii2h5part.cpp

clean:
	rm tomascii2h5part 2>&1 >/dev/null;\
	rm no_density   2>&1 >/dev/null | true


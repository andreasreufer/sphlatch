all: h5partDisplace h5partReadAttr h5partWriteAttr h5partCombineSolid h5partSlice h5partCombineA_ h5partCombine_E h5partCombine__

h5partDisplace:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -I../../src/old -lhdf5 -lz -lboost_program_options -o h5part_displace h5part_displace.cpp

h5partCombine__:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
	  -o h5part_combine__ h5part_combine.cpp

h5partCombineA_:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
	  -DSPHLATCH_ANEOS \
	  -o h5part_combineA_ h5part_combine.cpp

h5partCombine_E:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
	  -DSPHLATCH_GRAVITY_EPSSMOOTHING \
	  -o h5part_combine_E h5part_combine.cpp

h5partCombineSolid:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -DSPHLATCH_SOLID -o h5part_combine_SOLID h5part_combine.cpp

h5partReadAttr:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_readattr h5part_readattr.cpp

h5partWriteAttr:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src  -I../../src/old -lhdf5 -lz -lboost_program_options -o h5part_writeattr h5part_writeattr.cpp

h5partSlice:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_slice h5part_slice.cpp

icePeakP:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o ice_peakp ice_peakp.cpp

cutSphere:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -I../../src \
	  -fopenmp \
	  -o cut_sphere cut_sphere.cpp


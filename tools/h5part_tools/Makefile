all: h5partDisplace h5partCombine h5partReadAttr h5partWriteAttr h5partCombineSolid h5partSlice

h5partDisplace:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -I../../src/old -lhdf5 -lz -lboost_program_options -o h5part_displace h5part_displace.cpp

h5partCombine:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_combine h5part_combine.cpp

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

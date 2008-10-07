all: h5partDisplace h5partCombine h5partReadAttr h5partWriteAttr

h5partDisplace:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_displace h5part_displace.cpp

h5partCombine:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_combine h5part_combine.cpp

h5partReadAttr:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_readattr h5part_readattr.cpp

h5partWriteAttr:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_writeattr h5part_writeattr.cpp


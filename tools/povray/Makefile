all: hdf5topovrayCut
hdf5topovray:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o hdf52povray hdf52povray.cpp
	
hdf5topovrayCut:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o hdf52povray_cut hdf52povray_cut.cpp

ascii2povray:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lboost_program_options -o ascii2povray ascii2povray.cpp


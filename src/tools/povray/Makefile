all: hdf5topovrayCut
hdf5topovray:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o hdf52povray hdf52povray.cpp
	
hdf5topovrayCut:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o hdf52povray_cut hdf52povray_cut.cpp

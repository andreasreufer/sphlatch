all: tomsasciitoh5part hdf5toh5part
cdattoh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o cdat2h5part cdat2h5part.cpp
	
tomsasciitoh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o tomascii2h5part tomascii2h5part.cpp

hdf5toh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I./h5part -lhdf5 -lz -lboost_program_options -o hdf52h5part hdf52h5part.cpp

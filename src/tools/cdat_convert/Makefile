all: cdat2ascii cdat2float cdat2hdf5 cdat2vtk tipsyascii2cdat tomascii2cdat
cdat2ascii:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lboost_program_options -o cdat2ascii cdat2ascii.cpp
cdat2float:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lboost_program_options -o cdat2float cdat2float.cpp
cdat2hdf5:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lboost_program_options -o cdat2hdf5 cdat2hdf5.cpp
cdat2vtk:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lboost_program_options -o cdat2vtk cdat2vtk.cpp
#hdf52cdat:
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lboost_program_options -o hdf52cdat hdf52cdat.cpp

tipsyascii2cdat:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lboost_program_options -o tipsyascii2cdat tipsyascii2cdat.cpp

tomascii2cdat:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lboost_program_options -o tomascii2cdat tomascii2cdat.cpp


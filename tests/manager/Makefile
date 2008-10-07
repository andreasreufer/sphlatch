all: particle_manager communication_manager io_manager

particle_manager:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lboost_program_options -lhdf5 -lz -o particle_manager_test particle_manager_test.cpp

communication_manager:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lboost_program_options -lhdf5 -lz -o communication_manager_test communication_manager_test.cpp

attributes:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -o attributes_test attributes_test.cpp

io_manager:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lboost_program_options -lhdf5 -lz -o io_manager_test io_manager_test.cpp


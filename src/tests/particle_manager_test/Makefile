all: particle_manager communication_manager

particle_manager:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lboost_program_options -lhdf5 -lz -o particle_manager_test particle_manager_test.cpp

communication_manager:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o communication_manager_test communication_manager_test.cpp


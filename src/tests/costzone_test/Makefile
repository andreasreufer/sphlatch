all: costzone_test
	
costzone_test:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options -o costzone_test costzone_test.cpp

clean: 
	touch costzone_test
	rm costzone_test 2>&1 >/dev/null

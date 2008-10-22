all: aneosTest tillotsonTest

aneosTest: aneos 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -o aneos_test eos_test.cpp

storePress: aneos 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -lgfortran libaneos.o \
				      -o store_press_T store_press.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -o store_press_I store_press.cpp

tillotsonTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_TILLOTSON \
				      -o tillotson_test eos_test.cpp

findrhoTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -o findrho_test findrho_test.cpp


aneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

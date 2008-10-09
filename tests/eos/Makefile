all: aneosTest tillotsonTest

aneosTest: aneos 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -o aneos_test eos_test.cpp

tillotsonTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_TILLOTSON \
				      -o tillotson_test eos_test.cpp

aneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

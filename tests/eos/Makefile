all: aneosTest tillotsonTest storePress

aneosTestA_: 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -o aneos_testA_ eos_test.cpp
aneosTestAT: 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -DSPHLATCH_LOGGER \
				      -lgfortran libaneos.o \
				      -o aneos_testAT eos_test.cpp
aneosTestM_: 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_MANEOS \
				      -lgfortran libmaneos.o \
				      -o aneos_testM_ eos_test.cpp
aneosTestMT: 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_MANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -DSPHLATCH_LOGGER \
				      -lgfortran libmaneos.o \
				      -o aneos_testMT eos_test.cpp

storePress: aneos 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -lgfortran libaneos.o \
				      -o store_press_T store_press.cpp

idealgasTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_TIMEDEP_ENTROPY \
				      -o idealgas_test idealgas_test.cpp
tillotsonTest:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_TILLOTSON \
				      -o tillotson_test eos_test.cpp

findrhoTest: aneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -lgfortran libaneos.o \
				      -o findrho_test findrho_test.cpp

partsTest_T: aneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -lgfortran libaneos.o \
				      -o parts_test_T eos_parts.cpp

partsTestAT:  aneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_ANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -lgfortran libaneos.o \
				      -o parts_testAT eos_parts.cpp

partsTestMT: maneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz \
				      -DSPHLATCH_MANEOS \
				      -DSPHLATCH_ANEOS_TABLE \
				      -lgfortran libmaneos.o \
				      -o parts_testMT eos_parts.cpp

maneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libmaneos.f -o libmaneos.o

aneos:
	$(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

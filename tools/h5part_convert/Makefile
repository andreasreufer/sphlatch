all: cdattoh5part tomsasciitoh5part h5partAddscalpos fuchsbintoh5part h5parttoparasph h5parttoparasphpeakp

cdattoh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o cdat2h5part cdat2h5part.cpp
	
cdattoascii:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lboost_program_options -o cdat2ascii cdat2ascii.cpp

tomsasciitoh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o tomascii2h5part tomascii2h5part.cpp

h5partAddscalpos:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o h5part_addscalpos h5part_addscalpos.cpp

fuchsbintoh5part:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options -o fuchsbin2h5part fuchsbin2h5part.cpp

h5parttoparasph: xdrlib
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -lhdf5 -lz -lboost_program_options xdrlib.o -o h5part2parasph h5part2parasph.cpp

h5parttoparasphpeakp: xdrlib
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -I../../src -DPEAKP -lhdf5 -lz -lboost_program_options xdrlib.o -o h5part2parasph_peakp h5part2parasph.cpp

xdrlib:
	$(CC) $(CFLAGS) -c ../../aux/libxdr/xdrlib.c -o xdrlib.o


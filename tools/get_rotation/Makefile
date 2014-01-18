all: getRotation

getRotation:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -fopenmp \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_ANEOS_TRACKORIG \
	   -o get_rotation get_rotation.cpp

clean:
	rm get_rotation?? 2>&1 >/dev/null; true


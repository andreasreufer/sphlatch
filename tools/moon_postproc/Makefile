all: moonPostProc

moonPostProc:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -I../../src \
				      -DSPHLATCH_ANEOS \
				      -o moon_pp moon_postproc.cpp

clean:
	rm moon_pp 2>&1 >/dev/null;\
	true

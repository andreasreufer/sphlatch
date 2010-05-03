all: loadbalance

loadbalance:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) \
	  -lhdf5 -I../../src \
	  -DSPHLATCH_GRAVITY \
	  -DSPHLATCH_GRAVITY \
	  -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	  -DSPHLATCH_TIMEDEP_SMOOTHING \
	  -DSPHLATCH_TIMEDEP_ENERGY \
	  -fopenmp \
	  -o loadBalance load_balance.cpp

clean:
	rm loadBalance 2>&1 >/dev/null; true


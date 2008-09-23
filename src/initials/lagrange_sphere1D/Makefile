all: lagrangeSolverTillotson lagrangeSolverANEOS
	
libaneos:
	$(FC) $(FFLAGS) -c libaneos.f

lagrangeSolverANEOS: libaneos
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_ANEOS -lgfortran libaneos.o \
				      -o lagrange_sphere1D_ANEOS \
				      lagrange_sphere1D.cpp
	
lagrangeSolverTillotson:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -lhdf5 -lz -lboost_program_options \
				      -DSPHLATCH_TILLOTSON \
				      -o lagrange_sphere1D_tillotson \
				      lagrange_sphere1D.cpp
	

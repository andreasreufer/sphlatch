all: simple_sph_GCSHUmD_ simple_sph_GCSHSm_F simple_sph_GMSHUmD_ simple_sph_GMSHSm_F

# letters:
#
# 1.    G  self-gravity
#       _  none
#
# 2./3. CS standard-SPH, density sum
#       CI  "          , integrate density
#       MS miscible SPH, density sum
#       MI  "          , integrate density
#
# 4.    H  integrate smoothing length
#       _  don't "
#
# 5.    U  integrate internal energy
#       S  integrate entropy
#       _  none
#
# 6.    m  M-ANEOS with tables
#       M  M-ANEOS w/o  tables
#       a  ANEOS with tables
#       A  ANEOS w/o tables
#       T  Tillotson EOS
#       _  ideal gas only
#
# 7.    C  find clumps
#       E  find clumps, remove escapees,
#       D  find clumps and disks, remove escapees
#       _  none
#
# 8.    F  friction
#       Z  no acceleration in Z-direction
#       _  none
#
	   

simple_sph_GCSHU___:
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -fopenmp \
	   -o simple_sph_GHU___ simple_sph.cpp

simple_sph_GCSHUa__: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -fopenmp \
	   -o simple_sph_GCSHUa__ simple_sph.cpp

simple_sph_GCSHUaC_: libaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS -lgfortran libaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_FIND_CLUMPS \
	   -fopenmp \
	   -o simple_sph_GCSHUaC_ simple_sph.cpp

simple_sph_GCSHUmD_: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_FIND_CLUMPS \
	   -DSPHLATCH_ESCAPEES \
	   -DSPHLATCH_LRDISK \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_TRACK_ACCP \
	   -fopenmp \
	   -o simple_sph_GCSHUmD_ simple_sph.cpp

simple_sph_GCSHSm_F: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENTROPY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_TRACK_ACCP \
	   -DSPHLATCH_FRICTION \
	   -fopenmp \
	   -o simple_sph_GCSHSm_F simple_sph.cpp

simple_sph_GMSHUmD_: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_MISCIBLE \
	   -DSPHLATCH_FIND_CLUMPS \
	   -DSPHLATCH_ESCAPEES \
	   -DSPHLATCH_LRDISK \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_TRACK_ACCP \
	   -fopenmp \
	   -o simple_sph_GMSHUmD_ simple_sph.cpp

simple_sph_GMSHSm_F: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_MISCIBLE \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENTROPY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_TRACK_ACCP \
	   -DSPHLATCH_FRICTION \
	   -fopenmp \
	   -o simple_sph_GMSHSm_F simple_sph.cpp

simple_sph_GMIHUmD_: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_GRAVITY \
	   -DSPHLATCH_GRAVITY_POTENTIAL \
	   -DSPHLATCH_GRAVITY_SPLINESMOOTHING \
	   -DSPHLATCH_MISCIBLE \
	   -DSPHLATCH_INTEGRATERHO \
	   -DSPHLATCH_FIND_CLUMPS \
	   -DSPHLATCH_ESCAPEES \
	   -DSPHLATCH_LRDISK \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_TRACK_ACCP \
	   -fopenmp \
	   -o simple_sph_GMIHUmD_ simple_sph.cpp

simple_sph__CSHUm_Z: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_TRACK_ACCP \
	   -DSPHLATCH_ZONLY  \
	   -fopenmp \
	   -o simple_sph__CSHUm_Z simple_sph.cpp

simple_sph__MSHUm_Z: libmaneos
	 $(CXX) $(CXXFLAGS) $(LDFLAGS) \
	   -lhdf5 -lz -I../../src \
	   -DSPHLATCH_MISCIBLE \
	   -DSPHLATCH_TIMEDEP_SMOOTHING \
	   -DSPHLATCH_TIMEDEP_ENERGY \
	   -DSPHLATCH_ANEOS \
	   -DSPHLATCH_MANEOS -lgfortran libmaneos.o \
	   -DSPHLATCH_ANEOS_TABLE \
	   -DSPHLATCH_TRACK_TMAX \
	   -DSPHLATCH_TRACK_PMAX \
	   -DSPHLATCH_TRACK_UAV \
	   -DSPHLATCH_ZONLY  \
	   -fopenmp \
	   -o simple_sph__MSHUm_Z simple_sph.cpp


libaneos:
	    $(FC) $(FFLAGS) -c ../../aux/libaneos/libaneos.f -o libaneos.o

libmaneos:
	    $(FC) $(FFLAGS) -c ../../aux/libaneos/libmaneos.f -o libmaneos.o


clean:
	rm simple_sph_* 2>&1 >/dev/null; true


#include <iostream>
#include <iomanip>

#define SPHLATCH_HDF5
#define SPHLATCH_ANEOS_TABLE  
#define SPHLATCH_LOGGER 

#include "typedefs.h"

#include "sph_fluid_particle.h"
class particle :
  public sphlatch::SPHfluidPart,
  public sphlatch::energyPart,
  public sphlatch::ANEOSPart
{};

typedef particle partT;

#include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT> eosT;

int main(int argc, char* argv[])
{
  eosT& EOS(eosT::instance());  

#ifdef SPHLATCH_ANEOS
  EOS.storeTableS("aneos_tables.hdf5",2);
  EOS.storeTableU("aneos_tables.hdf5",2);
  EOS.storeTableS("aneos_tables.hdf5",4);
  EOS.storeTableU("aneos_tables.hdf5",4);
  EOS.storeTableS("aneos_tables.hdf5",5);
  EOS.storeTableU("aneos_tables.hdf5",5);
#endif

#ifdef SPHLATCH_MANEOS
  EOS.storeTableS("aneos_tables.hdf5",1);
  EOS.storeTableU("aneos_tables.hdf5",1);
#endif
  
  return EXIT_SUCCESS;
}

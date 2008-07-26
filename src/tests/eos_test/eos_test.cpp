// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_PARALLEL

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "eos_idealgas.h"
typedef sphlatch::EOS<sphlatch::IdealGas> eos_idealgas_type;

#include "eos_tillotson.h"
typedef sphlatch::EOS<sphlatch::Tillotson> eos_tillotson_type;

#include "eos_mie_grueneisen.h"
typedef sphlatch::EOS<sphlatch::MieGrueneisen> eos_miegrueneisen_type;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  part_type& PartManager(part_type::instance());

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
  
  PartManager.setNoParts(1);
  PartManager.resizeAll();

  PartManager.rho(0) = 1.;
  PartManager.u(0) = 1.;

  PartManager.attributes["gamma"] = 1.4;
  PartManager.attributes["KperU"] = 2.068e4;

  eos_idealgas_type&  myEOS1(eos_idealgas_type::instance());
  eos_idealgas_type&  myEOS2(eos_idealgas_type::instance());
  eos_tillotson_type& myEOS3(eos_tillotson_type::instance());

  const size_t i = 0;
  std::cout << myEOS1.getPressure(i) << "\n";
  std::cout << myEOS2.getPressure(i) << "\n";
  std::cout << myEOS3.getPressure(i) << "\n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

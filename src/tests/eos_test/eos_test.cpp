// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_PARALLEL

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

//#include "eos_idealgas.h"
//typedef sphlatch::EOS<sphlatch::IdealGas> eos_idealgas_type;

#include "eos_tillotson.h"
//typedef sphlatch::EOS<sphlatch::Tillotson> eos_tillotson_type;
typedef sphlatch::Tillotson eos_tillotson_type;

//#include "eos_mie_grueneisen.h"
//typedef sphlatch::EOS<sphlatch::MieGrueneisen> eos_miegrueneisen_type;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  part_type& PartManager(part_type::instance());

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
  
  PartManager.setNoParts(1);
  PartManager.resizeAll();


  PartManager.attributes["gamma"] = 1.4;
  PartManager.attributes["KperU"] = 2.068e4;

  //eos_idealgas_type&  myEOS1(eos_idealgas_type::instance());
  //eos_idealgas_type&  myEOS2(eos_idealgas_type::instance());
  eos_tillotson_type& myEOS3(eos_tillotson_type::instance());
  eos_tillotson_type& myEOS2(eos_tillotson_type::instance());
  
  //sphlatch::Tillotson myEOS3;

  
  sphlatch::valvectRefType rho(PartManager.rho);
  sphlatch::valvectRefType   u(PartManager.u);
  sphlatch::idvectRefType  mat(PartManager.mat);

  mat(0) = 1;
  
  const size_t idx = 0;
  rho(0) = 2.680e+00;
  u(0)   = 1.600e+11;
  for (size_t i = 0; i < 1024; i++)
  {
    std::cout << rho(0) << " " << u(0) << " " << myEOS3.getPressure(idx) << "\n";
    rho(0) -= 0.001;
    //u(0)   += 0.1e+11;
  }

  MPI::Finalize();
  return EXIT_SUCCESS;
}

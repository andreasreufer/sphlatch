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

#include "parasph_eos.cc"

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
  
  eosMat m;
  m.rho0  = 2.680e+00;
  m.A     = 1.800e+11;
  m.B     = 1.800e+11;
  m.a     = 0.5;
  m.b     = 1.3;
  m.alpha = 5.0;
  m.beta  = 5.0;
  m.u0    = 1.600e+11;
  m.Eiv   = 3.500e+10;
  m.Ecv   = 1.800e+11;

  const size_t idx = 0;
  rho(0) = 1.500e+00;
  //u(0)   = 3.000e+10;
  u(0)   = 0.000e+10;

  double rhom1, cs, T, press;
  int material = 0;

  sphlatch::valueType sphlatchP, sphlatchCs;

  for (size_t i = 0; i < 1024; i++)
  {
    rhom1 = 1. / rho(0);
    parasphTillotson( rho(0), rhom1, u(0), material, press, cs,
                      T, m );

    myEOS3.getPressCs(idx, sphlatchP, sphlatchCs);


    std::cout << rho(0) << " " << u(0) << " "
              << sphlatchP << " "
              << press << " "
              << sphlatchCs << " "
              << cs << "\n";
    
    //rho(0) += 0.005;
    u(0)   += 0.02e+10;
  }

  MPI::Finalize();
  return EXIT_SUCCESS;
}

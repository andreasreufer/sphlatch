// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_PARALLEL
#define SPHLATCH_ANEOS_TABLE

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

//#include "eos_idealgas.h"
//typedef sphlatch::IdealGas eos_type;

#ifdef SPHLATCH_TILLOTSON
#include "eos_tillotson.h"
#define SPHLATCH_EOS_DEFINED
typedef sphlatch::Tillotson eos_type;
#endif

#ifdef SPHLATCH_ANEOS
#include "eos_aneos.h"
#define SPHLATCH_EOS_DEFINED
typedef sphlatch::ANEOS eos_type;
#endif

//#include "eos_mie_grueneisen.h"
//typedef sphlatch::MieGrueneisen eos_type;

#include "parasph_eos.cc"

using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
#ifdef SPHLATCH_ANEOS
  PartManager.usePhase();
  PartManager.useTemperature();
#endif

  sphlatch::matrixRefType pos(PartManager.pos);
  sphlatch::valvectRefType rho(PartManager.rho);
  sphlatch::valvectRefType u(PartManager.u);
  sphlatch::idvectRefType mat(PartManager.mat);
  sphlatch::idvectRefType phase(PartManager.phase);

  IOManager.loadDump("initials.h5part");

  eos_type& EOS(eos_type::instance());

  /*eosMat dunite;
  dunite.rho0 = 3.320e+00;
  dunite.A = 1.290e+12;
  dunite.B = 1.290e+12;
  dunite.a = 0.5;
  dunite.b = 1.5;
  dunite.alpha = 5.0;
  dunite.beta = 5.0;
  dunite.u0 = 4.870e+12;
  dunite.Eiv = 4.720e+10;
  dunite.Ecv = 1.820e+11;

  eosMat iron;
  iron.rho0 = 7.860e+00;
  iron.A = 1.280e+12;
  iron.B = 1.050e+12;
  iron.a = 0.5;
  iron.b = 1.5;
  iron.alpha = 5.0;
  iron.beta = 5.0;
  iron.u0 = 9.500e+10;
  iron.Eiv = 1.420e+10;
  iron.Ecv = 8.450e+10;

  sphlatch::valueType sP = 0., sCs = 0.;

  double pP = 0., pCs = 0., rhom1 = 0., T = 0.;*/

  sphlatch::valueType sP = 0., sCs = 0.;
  const size_t noParts = PartManager.getNoLocalParts();
  //for (size_t i = 0; i < 1024; i++)
  for (size_t i = 0; i < noParts; i++)
    {
      /*rhom1 = 1. / rho(i);

      if (mat(i) == 4)
        parasphTillotson(rho(i), rhom1, u(i), mat(i),
                         pP, pCs, T, dunite);
      else if (mat(i) == 5)
        parasphTillotson(rho(i), rhom1, u(i), mat(i),
                         pP, pCs, T, iron);

      EOS(i, sP, sCs);

      std::cout << i << "   "
                << pos(i, X) << "   "
                << pos(i, Y) << "   "
                << pos(i, Z) << "   "
                << mat(i) << "   " //  5
                << u(i) << "   "
                << rho(i) << "   "
                << pP << "   "     //  8
                << pCs << "   "    //  9
                << sP << "   "     // 10
                << sCs << "   "    // 11
                << pP - sP << "    " // 12
                << pCs - sCs << "  " // 13
                << phase(i) << "\n"; // 14*/
     
      mat(i) = 4;
      
      rho(i) = 8.70359;
      u(i) = 0.403702;

      EOS(i, sP, sCs);

      std::cout << mat(i) << "\t"
                << sP << "\t"
                << sCs << "\n";
    }

  MPI::Finalize();
  return EXIT_SUCCESS;
}

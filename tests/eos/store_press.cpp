// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_PARALLEL
//#define SPHLATCH_ANEOS_TABLE

#define SPHLATCH_LOGGER

#include <iostream>
#include <iomanip>
#include <boost/assign/std/vector.hpp>

#include "typedefs.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "log_manager.h"
typedef sphlatch::LogManager log_type;

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

using namespace sphlatch::vectindices;
using namespace boost::assign;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());
  log_type& Logger(log_type::instance());

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
#ifdef SPHLATCH_ANEOS
  PartManager.usePhase();
  PartManager.useTemperature();
#endif

  sphlatch::matrixRefType pos(PartManager.pos);
  
  sphlatch::valvectRefType rho(PartManager.rho);
  sphlatch::valvectRefType   u(PartManager.u);
  sphlatch::valvectRefType  cs(PartManager.cs);
  sphlatch::valvectRefType   p(PartManager.p);
#ifdef SPHLATCH_ANEOS
  sphlatch::valvectRefType   T(PartManager.T);
#endif

  sphlatch::idvectRefType    id(PartManager.id);
  sphlatch::idvectRefType   mat(PartManager.mat);
#ifdef SPHLATCH_ANEOS
  sphlatch::idvectRefType phase(PartManager.phase);
#endif

  IOManager.loadDump("input.h5part");

  PartManager.attributes["rhomin"] = 1.e-6;
  PartManager.attributes["rhomax"] = 15.;

  PartManager.attributes["umin"] = 1.e6;
  PartManager.attributes["umax"] = 1.e14;

  eos_type& EOS(eos_type::instance());

  ///
  /// first "slow" run
  ///
  const size_t noParts = PartManager.getNoLocalParts();
  for (size_t i = 0; i < noParts; i++)
    {
      EOS(i, p(i), cs(i));
    }

  Logger << "first EOS run (slow)";
  
  for (size_t i = 0; i < noParts; i++)
    {
      EOS(i, p(i), cs(i));
    }
  
  Logger << "second EOS run (fast)";

  sphlatch::quantsType saveQuants;
  saveQuants.ints += &id, &mat;
#ifdef SPHLATCH_ANEOS
  saveQuants.ints += &phase;
#endif
  saveQuants.vects += &pos;
  saveQuants.scalars += &p, &cs, &rho, &u;
#ifdef SPHLATCH_ANEOS
  saveQuants.scalars += &T;
#endif

  IOManager.saveDump("output.h5part", saveQuants);

  MPI::Finalize();
  return EXIT_SUCCESS;
}

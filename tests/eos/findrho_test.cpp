// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_PARALLEL

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

//#include "eos_idealgas.h"
//typedef sphlatch::IdealGas eos_type;

#include "eos_tillotson.h"
typedef sphlatch::Tillotson eos_type;

using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  eos_type& EOS(eos_type::instance());

  std::cout << EOS.findRho( 1.72e9, 4, 0., 1., 1.0, 10. ) << "\n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

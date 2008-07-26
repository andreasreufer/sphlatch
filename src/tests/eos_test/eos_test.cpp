// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#define SPHLATCH_PARALLEL

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "eos_idealgas.h"
#include "eos_tillotson.h"
typedef sphlatch::EOS<sphlatch::IdealGas> eos1_type;
typedef sphlatch::EOS<sphlatch::IdealGas> eos2_type;
typedef sphlatch::EOS<sphlatch::Tillotson> eos3_type;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  eos1_type& myEOS1(eos1_type::instance());
  eos2_type& myEOS2(eos2_type::instance());
  eos3_type& myEOS3(eos3_type::instance());

  myEOS1.pressure();
  myEOS2.pressure();
  myEOS3.pressure();

  MPI::Finalize();
  return EXIT_SUCCESS;
}

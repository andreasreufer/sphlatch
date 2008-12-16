// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define GRAVITY
#define TREEORDER
//#define RSORDER

//#define SPHLATCH_RANKSPACESERIALIZE

//#define BFCHECK
//#define CHECK_TREE
//#define CHECK_RANKSPACE

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <boost/assign/std/vector.hpp>

#include "typedefs.h"
typedef sphlatch::fType fType;

#include "mpi.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager com_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

using namespace boost::assign;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  part_type& PartManager(part_type::instance());

  fType& time(PartManager.attributes["time"]);

  std::cout << time << "\n";
  time = 12.;
  std::cout << PartManager.attributes["time"] << "\n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}



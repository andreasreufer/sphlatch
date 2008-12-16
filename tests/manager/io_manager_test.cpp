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

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include <boost/assign/std/vector.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::fType fType;
typedef sphlatch::stringListType stringListType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::matrixType matrixType;

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
  double stepStartTime, logStartTime = MPI_Wtime();


  po::options_description Options("Global Options");
  Options.add_options() ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file");

  po::positional_options_description POD;
  POD.add("input-file", 1);

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(POD).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  if (!VMap.count("input-file") or
      !VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  com_type& ComManager(com_type::instance());
  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  const size_t myDomain = ComManager.getMyDomain();
  const size_t noDomains = ComManager.getNoDomains();

  std::string inputFileName = VMap[ "input-file"].as<std::string>();
  std::string outputFileName = VMap["output-file"].as<std::string>();

  stringListType vars = IOManager.discoverVars(inputFileName);

  stringListType::iterator varItr = vars.begin();
  stringListType::iterator varEnd = vars.end();

  while( varItr != varEnd )
  {
    std::cout << "found " << *varItr << "\n";
    varItr++;
  }

  PartManager.useBasicSPH();

  IOManager.loadDump(inputFileName);

  sphlatch::quantsType saveQuants;
  saveQuants.vects += &PartManager.pos;
  saveQuants.scalars += &PartManager.m;

  IOManager.saveDump(outputFileName, saveQuants);
  valvectType vect(10);
  matrixType matr(10, 3);
  
  for (size_t i = 0; i < vect.size(); i++)
  {
    vect(i) = i;
    matr(i, 0) = i + 11;
    matr(i, 1) = i + 21;
    matr(i, 2) = i + 31;
  }
  IOManager.savePrimitive( vect, "haba", "test.hdf5" );
  IOManager.savePrimitive( matr, "haba2", "test.hdf5" );

  MPI::Finalize();
  return EXIT_SUCCESS;
}



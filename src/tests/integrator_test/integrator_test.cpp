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
#define RSORDER

//#define SPHLATCH_RANKSPACESERIALIZE

//#define BFCHECK
//#define CHECK_TREE
//#define CHECK_RANKSPACE

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include <boost/assign/std/vector.hpp>

#include <boost/mpl/vector_c.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
//typedef sphlatch::valueType valueType;
//typedef sphlatch::partsIndexVectType partsIndexVectType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager com_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include "integrator_verlet.h"

#include <boost/progress.hpp>
#include <vector>

// tree stuff
#include "bhtree.h"

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

  using namespace boost::assign;
  using namespace sphlatch;

  io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  com_type& ComManager(com_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  const size_t myDomain = ComManager.getMyDomain();

  size_t noParts, noGhosts, noTotParts;

  std::string inputFileName = VMap[ "input-file"].as<std::string>();
  std::string outputFileName = VMap["output-file"].as<std::string>();

  PartManager.useGravity();
  
  IOManager.loadDump(inputFileName);

  // set up logging stuff
  std::string logFilename = outputFileName + "_rank000";
  std::string rankString = boost::lexical_cast<std::string>(myDomain);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "    start log\n";

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);

  idvectRefType id(PartManager.id);

  quantsType exchQuants;
  exchQuants.vects += &pos, &vel;
  exchQuants.ints += &id;
  exchQuants.scalars += &m;

  stepStartTime = MPI_Wtime();
  /// domain decomposition and exchange of particles
  CostZone.createDomainPartsIndex();
  ComManager.exchange(CostZone.domainPartsIndex,
                      CostZone.getNoGhosts(),
                      exchQuants);

  /// prepare ghost stuff
  ComManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  ComManager.sendGhosts(exchQuants);

  noParts = PartManager.getNoLocalParts();
  noGhosts = PartManager.getNoGhostParts();
  noTotParts = PartManager.getNoTotalParts();

  size_t curIndex = 0;
  double treesearchTime;

  MetaIntegrator<Verlet> myIntegrator;

  myIntegrator.regIntegration(m, m, m);


  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}



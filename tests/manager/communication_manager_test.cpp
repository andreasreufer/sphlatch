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

//#define SPHLATCH_RANKSPACESERIALIZE

// enable checking of bounds for the neighbour lists
#define SPHLATCH_CHECKNONEIGHBOURS

// enable selfgravity?
//#define SPHLATCH_GRAVITY

// time dependent energy?
//#define SPHLATCH_TIMEDEP_ENERGY

// time dependent smoothing length?
//#define SPHLATCH_TIMEDEP_SMOOTHING

// shocktube boundary conditions?
//#define SPHLATCH_SHOCKTUBE

// integrate rho instead of using the SPH sum?
//#define SPHLATCH_INTEGRATERHO

// linear velocity damping term?
//#define SPHLATCH_FRICTION

// do we need the velocity divergence?
#ifdef SPHLATCH_TIMEDEP_ENERGY
#ifndef SPHLATCH_VELDIV
#define SPHLATCH_VELDIV
#endif
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
#ifndef SPHLATCH_VELDIV
#define SPHLATCH_VELDIV
#endif
#endif

#ifdef SPHLATCH_INTEGRATERHO
#ifndef SPHLATCH_VELDIV
#define SPHLATCH_VELDIV
#endif
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

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
typedef sphlatch::fType fType;
typedef sphlatch::fRefType fRefType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::valvectRefType valvectRefType;

typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::partsIndexVectType partsIndexVectType;
typedef sphlatch::quantsType quantsType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager comm_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include "log_manager.h"
typedef sphlatch::LogManager log_type;

#include <vector>

using namespace sphlatch::vectindices;
using namespace boost::assign;

int main(int argc, char* argv[])
{
  ///
  /// parse program options
  ///
  po::options_description Options("<input-file>\n ... or use options");

  Options.add_options()
  ("input-file,i", po::value<std::string>(), "input file");

  po::positional_options_description posDesc;
  posDesc.add("input-file", 1);

  po::variables_map poMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(posDesc).run(), poMap);
  po::notify(poMap);

  /*if (!poMap.count("input-file"))
    {
      std::cerr << Options << "\n";
      return EXIT_SUCCESS;
    }*/

  ///
  /// everythings set, now start the parallel envirnoment
  ///
  MPI::Init(argc, argv);

  ///
  /// instantate managers
  ///
  //io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  log_type& Logger(log_type::instance());

  ///
  /// some simulation parameters
  /// attributes will be overwritten, when defined in file
  ///
  //std::string loadDumpFile = poMap["input-file"].as<std::string>();

  ///
  /// define what we're doing
  ///
  PartManager.useBasicSPH();

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  valvectRefType m(PartManager.m);
  idvectRefType id(PartManager.id);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos;
  CommManager.exchangeQuants.scalars += &m;
  CommManager.exchangeQuants.ints += &id;

  ///
  /// log program compilation time
  ///
  Logger.stream << "executable compiled from " << __FILE__
                << " on " << __DATE__
                << " at " << __TIME__ << "\n\n"
                << "    features: \n"
                << "     basic SPH\n";
  Logger.flushStream();

  ///
  /// load particles
  ///
  //IOManager.loadDump(loadDumpFile);
  //Logger.stream << "loaded " << loadDumpFile;
  //Logger.flushStream();

  ///
  /// generate particles
  ///
  const size_t noGenParts = 30;
  PartManager.setNoParts(noGenParts);
  PartManager.resizeAll();
  for (size_t i = 0; i < noGenParts; i++)
    {
      id(i) = noGenParts * myDomain + i + 1;
      m(i) = noGenParts * myDomain + i + 1;

      pos(i, X) = static_cast<sphlatch::fType>(rand()) / RAND_MAX;
      pos(i, Y) = static_cast<sphlatch::fType>(rand()) / RAND_MAX;
      pos(i, Z) = static_cast<sphlatch::fType>(rand()) / RAND_MAX;

    }

  ///
  /// exchange particle data
  ///
  CostZone.createDomainPartsIndex();
  Logger << "new domain decomposition";
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());

  ///
  /// prepare ghost sends
  ///
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  Logger.stream << "distributed particles: "
                << PartManager.getNoLocalParts() << " parts. & "
                << PartManager.getNoGhostParts() << " ghosts";
  Logger.flushStream();

  const size_t noParts = PartManager.getNoLocalParts();
  //const size_t noTotParts = noParts + PartManager.getNoGhostParts();

  ///
  /// send ghosts to other domains
  ///
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);

  std::cout << myDomain << ": id " << id << "\n\n";

  CommManager.barrier();

  std::cout << myDomain << ": m " << m << "\n\n";
  
  CommManager.barrier();

  std::cout << myDomain << ": pos " << pos << "\n\n";

  ///
  /// now move the particle a little bit
  ///
  /*for (size_t i = 0; i < noParts; i++)
    {
      pos(i, X) += 0.05 * (
        static_cast<sphlatch::fType>(rand()) / RAND_MAX - 0.5);
      pos(i, Y) += 0.05 * (
        static_cast<sphlatch::fType>(rand()) / RAND_MAX - 0.5);
      pos(i, Z) += 0.05 * (
        static_cast<sphlatch::fType>(rand()) / RAND_MAX - 0.5);

      ///
      /// this ensures the particles to stay in a periodic boundary box
      ///
      if (pos(i, X) > 1.00)
        pos(i, X) -= 1.00;
      if (pos(i, Y) > 1.00)
        pos(i, Y) -= 1.00;
      if (pos(i, Z) > 1.00)
        pos(i, Z) -= 1.00;

      if (pos(i, X) < 0.00)
        pos(i, X) += 1.00;
      if (pos(i, Y) < 0.00)
        pos(i, Y) += 1.00;
      if (pos(i, Z) < 0.00)
        pos(i, Z) += 1.00;
    }*/

  Logger << "move particles by 0.05 at most";

  MPI::Finalize();
  return EXIT_SUCCESS;
}



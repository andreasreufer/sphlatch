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

// enable checking of bounds for the neighbour lists
//#define SPHLATCH_CHECKNONEIGHBOURS

// enable selfgravity?
//#define SPHLATCH_GRAVITY

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
typedef sphlatch::valueType valueType;
typedef sphlatch::valueRefType valueRefType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::zerovalvectType zerovalvectType;
typedef sphlatch::valvectRefType valvectRefType;

typedef sphlatch::particleRowType particleRowType;

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

#include <boost/progress.hpp>
#include <vector>

/// tree stuff
#include "bhtree.h"

using namespace sphlatch::vectindices;
using namespace boost::assign;



int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  ///
  /// instantate managers
  ///
  io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  log_type& Logger(log_type::instance());

  ///
  /// some simulation parameters
  /// attributes will be overwritten, when defined in file
  ///
  PartManager.attributes["gravconst"] = 1.0;
  PartManager.attributes["gravtheta"] = 0.7;

  ///
  /// define what we're doing and load data
  ///
  PartManager.useGravity();
  IOManager.loadDump("initial.h5part");

  ///
  /// some useful references
  ///
  idvectRefType id(PartManager.id);
  
  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType eps(PartManager.eps);

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &eps;
  CommManager.exchangeQuants.ints += &id;


  ///
  /// exchange particles
  ///
  CostZone.createDomainPartsIndex();
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
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();

  /// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(eps); // << eps is not used for interacting partners!
  Logger << " sent to ghosts: pos, vel, id, m, eps";

  const valueType gravTheta = PartManager.attributes["gravtheta"];
  const valueType gravConst = PartManager.attributes["gravconst"];

  //sphlatch::BHtree<sphlatch::Monopoles> Tree(gravTheta,
  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
                                               gravConst,
                                               CostZone.getDepth(),
                                               CostZone.getCenter(),
                                               CostZone.getSidelength()
                                               );

  ///
  /// fill up tree, determine ordering, calculate multipoles
  ///
  for (size_t i = 0; i < noTotParts; i++)
    {
      Tree.insertParticle(i);
    }

  const size_t myDomain = CommManager.getMyDomain();
  std::string domString = "domain";
  domString.append( boost::lexical_cast<std::string>( myDomain ) );
  std::string treeDumpfile;
  
  treeDumpfile = "treeDump_preMP.";
  treeDumpfile.append( domString );
  treeDumpfile.append( ".dot" );
  Tree.treeDOTDump( treeDumpfile );
  Logger.stream << "tree dumped to " << treeDumpfile;
  Logger.flushStream();
 
  ///
  /// determine particle order and calculate multipoles
  ///
  Tree.detParticleOrder();
  Tree.calcMultipoles();
  Logger << "Tree ready";

  ///
  /// inspect the tree
  ///
  std::string toptreeDumpfile = "toptreeDump.";
  toptreeDumpfile.append( domString );
  Tree.toptreeDump( toptreeDumpfile );
  Logger.stream << "toptree dumped to " << toptreeDumpfile;
  Logger.flushStream();
  
  treeDumpfile = "treeDump_postMP.";
  treeDumpfile.append( domString );
  treeDumpfile.append( ".dot" );
  Tree.treeDOTDump( treeDumpfile );
  Logger.stream << "tree dumped to " << treeDumpfile;
  Logger.flushStream();


  ///
  /// calc gravity
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      const size_t curIdx = Tree.particleOrder[i];
      Tree.calcGravity(curIdx);
    }
  Logger << "gravity calculated";

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &eps;
  saveQuants.ints += &id;

  IOManager.saveDump("saveDump.h5part", saveQuants);

  MPI::Finalize();
  return EXIT_SUCCESS;
}



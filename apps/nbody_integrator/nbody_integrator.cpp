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
//#define SPHLATCH_CHECKNONEIGHBOURS

// enable selfgravity?
#define SPHLATCH_GRAVITY

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

#include "integrator_verlet.h"

#include <boost/progress.hpp>
#include <vector>

/// tree stuff
#include "bhtree.h"

using namespace sphlatch::vectindices;
using namespace boost::assign;

valueType timeStep()
{
  part_type& PartManager(part_type::instance());
  //comm_type& CommManager(comm_type::instance());
  log_type& Logger(log_type::instance());
  io_type& IOManager(io_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType eps(PartManager.eps);

  idvectRefType id(PartManager.id);

  valueRefType time(PartManager.attributes["time"]);
  int& step(PartManager.step);
  //int& substep(PartManager.substep);
  //const size_t noParts = PartManager.getNoLocalParts();
  //const size_t myDomain = CommManager.getMyDomain();

  ///
  /// timestep dictated by acceleration
  ///
  const valueType dtA = 0.002;

  ///
  /// distance to next saveItrvl
  ///
  //const valueType saveItrvl = 0.1;
  const valueType saveItrvl = 0.002;
  std::string fileName = "saveDump000000.h5part";
  const valueType nextSaveTime = (floor((time / saveItrvl) + 1.e-6)
                                  + 1.) * saveItrvl;
  valueType dtSave = nextSaveTime - time;

  ///
  /// determine global minimum.
  /// by parallelly minimizing the timesteps, we
  /// can estimate which ones are dominant
  ///
  valueType dtGlob = std::numeric_limits<valueType>::max();
  
  dtGlob = dtA < dtGlob ? dtA : dtGlob;
  dtGlob = dtSave < dtGlob ? dtSave : dtGlob;

  Logger.stream << "dtA: " << dtA
                << " dtSave: " << dtSave
                << "   dtGlob: " << dtGlob;
  Logger.flushStream();

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos, &vel, &acc;
  saveQuants.scalars += &m, &eps;
  saveQuants.ints += &id;

  const valueType curSaveTime = (floor((time / saveItrvl) + 1.e-9))
                                * saveItrvl;
  if (fabs(curSaveTime - time) < 1.e-9)
    {
      Logger << "save dump";

      std::string stepString = boost::lexical_cast<std::string>(step);

      fileName.replace(fileName.size() - 7 - stepString.size(),
                       stepString.size(), stepString);

      IOManager.saveDump(fileName, saveQuants);
    }

  return dtGlob;
}

void derivate()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  //io_type& IOManager(io_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  log_type& Logger(log_type::instance());

  matrixRefType pos(PartManager.pos);

  valvectRefType m(PartManager.m);
  valvectRefType eps(PartManager.eps);

  idvectRefType id(PartManager.id);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  //int& step(PartManager.step);
  //int& substep(PartManager.substep);
  //const size_t myDomain = CommManager.getMyDomain();

/// little helper vector to zero a 3D quantity
  const zerovalvectType zero(3);

/// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);
  CommManager.sendGhosts(eps); // << eps is not used for interacting partners!
  Logger << " sent to ghosts: pos, id, m, eps";

  const valueType gravTheta = PartManager.attributes["gravtheta"];
  const valueType gravConst = PartManager.attributes["gravconst"];

  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
                                               gravConst,
                                               CostZone.getDepth(),
                                               CostZone.getCenter(),
                                               CostZone.getSidelength()
                                               );

#ifdef SPHLATCH_GRAVITY
  ///
  /// fill up tree, determine ordering, calculate multipoles
  ///
  for (size_t i = 0; i < noTotParts; i++)
    {
      Tree.insertParticle(i);
    }

  Tree.detParticleOrder();
  Tree.calcMultipoles();
  Logger << "Tree ready";

  for (size_t i = 0; i < noParts; i++)
    {
      const size_t curIdx = Tree.particleOrder[i];
      Tree.calcGravity(curIdx);
    }
  Logger << "gravity calculated";
#endif
};


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

  std::string loadDumpFile = "initial.h5part";
  const valueType maxTime = 0.01;

  ///
  /// define what we're doing
  ///
  PartManager.useGravity();

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);
  valvectRefType eps(PartManager.eps);

  idvectRefType id(PartManager.id);

  int& step(PartManager.step);
  valueRefType time(PartManager.attributes["time"]);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m, &eps;
  CommManager.exchangeQuants.ints += &id;
  
  ///
  /// instantate the MetaIntegrator
  ///
  sphlatch::VerletMetaIntegrator Integrator(derivate, timeStep);
  //sphlatch::PredCorrMetaIntegrator Integrator(derivate, timeStep);

  ///
  /// register spatial, energy and smoothing length integration
  ///
  Integrator.regIntegration(pos, vel, acc);

  IOManager.loadDump(loadDumpFile);
  Logger.stream << "loaded " << loadDumpFile;
  Logger.flushStream();

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

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap();
  Logger << "integrator bootstrapped";

  ///
  /// the integration loop
  ///
  //valueType nextSaveTime = time;
  while (time <= maxTime)
    {
      ///
      /// exchange particles and ghosts
      /// \todo: only do this, when necessary
      ///
      CostZone.createDomainPartsIndex();
      CommManager.exchange(CostZone.domainPartsIndex,
                           CostZone.getNoGhosts());

      CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
      Logger.stream << "distributed particles: "
                    << PartManager.getNoLocalParts() << " parts. & "
                    << PartManager.getNoGhostParts() << " ghosts";
      Logger.flushStream();

      ///
      /// integrate
      ///
      Integrator.integrate();

      Logger.stream << "finished step " << step << ", now at t = " << time;
      Logger.flushStream();
      Logger.zeroRelTime();

      if (myDomain == 0)
        {
          std::cout << "t = " << std::fixed << std::right
                    << std::setw(12) << std::setprecision(6)
                    << time << " (" << step << ")\n";
        }
    }

  MPI::Finalize();
  return EXIT_SUCCESS;
}



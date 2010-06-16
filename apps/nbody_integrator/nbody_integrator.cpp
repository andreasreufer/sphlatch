// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable Logger
#define SPHLATCH_LOGGER

//#define SPHLATCH_RANKSPACESERIALIZE

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
typedef sphlatch::fType fType;
typedef sphlatch::fRefType fRefType;

typedef sphlatch::valvectType valvectType;
typedef sphlatch::valvectRefType valvectRefType;

typedef sphlatch::countsType countsType;

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
typedef sphlatch::VerletMetaIntegrator integrator_type;

//#include "integrator_predcorr.h"
//typedef sphlatch::PredCorrMetaIntegrator integrator_type;


#include <boost/progress.hpp>
#include <vector>

#include "bhtree.h"
typedef sphlatch::BHtree<sphlatch::Quadrupoles> tree_type;

using namespace sphlatch::vectindices;
using namespace boost::assign;

fType timeStep()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  log_type& Logger(log_type::instance());
  io_type& IOManager(io_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);

#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif
  
  idvectRefType id(PartManager.id);

  fRefType time(PartManager.attributes["time"]);
  fRefType saveItrvl(IOManager.saveItrvl);

  int& step(PartManager.step);
  //int& substep(PartManager.substep);

  const size_t noParts = PartManager.getNoLocalParts();
  //const size_t myDomain = CommManager.getMyDomain();

  ///
  /// timestep criterion for acceleration
  /// ( see Wetzstein et. al 2008 )
  ///
  fType dtA = std::numeric_limits<fType>::max();
  /*for (size_t i = 0; i < noParts; i++)
    {
      const fType ai = sqrt(acc(i, X) * acc(i, X) +
                            acc(i, Y) * acc(i, Y) +
                            acc(i, Z) * acc(i, Z));

      if (ai > 0.)
        {
          const fType dtAi = sqrt(h(i) / ai);
          dtA = dtAi < dtA ? dtAi : dtA;
        }
    }
  dtA *= 0.5;*/
  dtA = 1.;
  CommManager.min(dtA);

  ///
  /// distance to next save time
  ///
  const fType dtSave = (floor((time / saveItrvl) + 1.e-6)
                        + 1.) * saveItrvl - time;

  ///
  /// determine global minimum.
  /// by parallelly minimizing the timesteps, we
  /// can estimate which ones are dominant
  ///
  fType dtGlob = std::numeric_limits<fType>::max();

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
  saveQuants.scalars += &m;
  saveQuants.ints += &id;

#ifdef SPHLATCH_GRAVITY
  saveQuants.scalars += &eps;
#endif
  const fType curSaveTime = (floor((time / saveItrvl) + 1.e-9))
                            * saveItrvl;
  if (fabs(curSaveTime - time) < 1.e-9)
    {
      std::string fileName = "dump";

      std::ostringstream stepSS;
      stepSS << step;

      ///
      /// pad step number to 7 numbers
      ///
      for (size_t i = stepSS.str().size(); i < 7; i++)
        {
          fileName.append("0");
        }
      fileName.append(stepSS.str());

      fileName.append("_T");
      std::ostringstream timeSS;
      timeSS << std::setprecision(4) << std::scientific << time;

      ///
      /// pad to a size of 11 characters
      /// (sign 1, mantissa 1, '.' 1, prec 3, 'e+' 2, exponent 3)
      ///
      for (size_t i = timeSS.str().size(); i < 11; i++)
        {
          fileName.append("0");
        }
      fileName.append(timeSS.str());
      fileName.append(".h5part");

      IOManager.saveDump(fileName, saveQuants);
      Logger.stream << "dump saved: " << fileName;
      Logger.flushStream();
    }

  return dtGlob;
}

void derivate()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  log_type& Logger(log_type::instance());

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  valvectRefType m(PartManager.m);

#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif

  idvectRefType id(PartManager.id);

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
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  int& step(PartManager.step);
  int& substep(PartManager.substep);
  //const size_t myDomain = CommManager.getMyDomain();

  /// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(m);
#ifdef SPHLATCH_GRAVITY
  CommManager.sendGhosts(eps); // << eps is not yet used for interacting partners!
#endif
  Logger << " sent to ghosts: pos, vel, id, m, h, eps";

  ///
  /// zero the derivatives
  ///
  for (size_t i = 0; i < noParts; i++)
    {
      acc(i, X) = 0.;
      acc(i, Y) = 0.;
      acc(i, Z) = 0.;
    }

#ifdef SPHLATCH_GRAVITY
  const fType gravTheta = PartManager.attributes["gravtheta"];
  const fType gravConst = PartManager.attributes["gravconst"];

  tree_type Tree(gravTheta, gravConst,
                 CostZone.getDepth(), CostZone.getCenter(),
                 CostZone.getSidelength());

  ///
  /// fill up tree, determine ordering, calculate multipoles
  ///
  for (size_t k = 0; k < noTotParts; k++)
    {
      Tree.insertParticle(k);
    }

  Tree.detParticleOrder();
  Tree.calcMultipoles();
  Logger << "Tree ready";

  for (size_t k = 0; k < noParts; k++)
    {
      const size_t i = Tree.particleOrder[k];
      Tree.calcAcc(i);
    }
  Logger << "gravity calculated";
#endif

};

int main(int argc, char* argv[])
{
  ///
  /// parse program options
  ///
  po::options_description Options(
    "<input-file> <save-time> <stop-time>\n ... or use options");

  Options.add_options()
  ("input-file,i", po::value<std::string>(),
   "input file")
  ("save-time,s", po::value<fType>(),
   "save dumps when (time) modulo (save time) = 0.")
  ("stop-time,S", po::value<fType>(),
   "stop simulaton at this time");

  po::positional_options_description posDesc;
  posDesc.add("input-file", 1);
  posDesc.add("save-time", 1);
  posDesc.add("stop-time", 1);

  po::variables_map poMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(posDesc).run(), poMap);
  po::notify(poMap);

  if (!poMap.count("input-file") ||
      !poMap.count("save-time") ||
      !poMap.count("stop-time"))
    {
      std::cerr << Options << "\n";
      return EXIT_FAILURE;
    }

  ///
  /// everythings set, now start the parallel envirnoment
  ///
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
  std::string loadDumpFile = poMap["input-file"].as<std::string>();
  const fType maxTime = poMap["stop-time"].as<fType>();
  IOManager.saveItrvl = poMap["save-time"].as<fType>();

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
#ifdef SPHLATCH_GRAVITY
  valvectRefType eps(PartManager.eps);
#endif

  idvectRefType id(PartManager.id);

  int& step(PartManager.step);
  fRefType time(PartManager.attributes["time"]);

  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &m;
  CommManager.exchangeQuants.ints += &id;
#ifdef SPHLATCH_GRAVITY
  CommManager.exchangeQuants.scalars += &eps;
#endif

  ///
  /// instantate the MetaIntegrator
  ///
  sphlatch::VerletMetaIntegrator Integrator(derivate, timeStep);
  //sphlatch::PredCorrMetaIntegrator Integrator(derivate, timeStep);

  ///
  /// register spatial, energy and smoothing length integration
  ///
  Integrator.regIntegration(pos, vel, acc);

  ///
  /// log program compilation time
  ///
  Logger.stream << "executable compiled from " << __FILE__
                << " on " << __DATE__
                << " at " << __TIME__ << "\n\n"
                << "    features: \n"
#ifdef SPHLATCH_GRAVITY
                << "     gravity  \n"
#endif
                << "";
  Logger.flushStream();

  ///
  /// load particles
  ///
  IOManager.loadDump(loadDumpFile);
  Logger.stream << "loaded " << loadDumpFile;
  Logger.flushStream();

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap();
  Logger << "integrator bootstrapped";

  ///
  /// the integration loop
  ///
  while (time <= maxTime)
    {
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



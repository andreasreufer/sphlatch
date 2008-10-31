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

// remove particles on escaping orbits?
//#define SPHLATCH_REMOVEESCAPING

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
typedef sphlatch::valueType valueType;
typedef sphlatch::valueRefType valueRefType;

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

#include "kernel_cubicspline.h"
typedef sphlatch::CubicSpline2D kernel_type;

#include "rankspace.h"
typedef sphlatch::Rankspace neighsearch_type;

using namespace sphlatch::vectindices;
using namespace boost::assign;



int main(int argc, char* argv[])
{
  ///
  /// parse program options
  ///
  po::options_description Options(
    "<input-file> \n ... or use options");

  Options.add_options()
  ("input-file,i", po::value<std::string>(),
    "input file");

  po::positional_options_description posDesc;
  posDesc.add("input-file", 1);

  po::variables_map poMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(posDesc).run(), poMap);
  po::notify(poMap);

  if (!poMap.count("input-file") )
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
  log_type& Logger(log_type::instance());
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  io_type& IOManager(io_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  ///
  /// define what we're doing
  ///
  PartManager.useBasicSPH();
  
  ///
  /// some refs
  ///
  matrixRefType pos(PartManager.pos);
  valvectRefType h(PartManager.h);
  valvectRefType rho(PartManager.rho);
  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos;
  CommManager.exchangeQuants.scalars += &h;
  CommManager.exchangeQuants.ints += &id;

  ///
  /// load particles
  ///
  std::string loadDumpFile = poMap["input-file"].as<std::string>();
  IOManager.loadDump(loadDumpFile);
  Logger.stream << "loaded " << loadDumpFile;
  Logger.flushStream();

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

  /// send ghosts to other domains
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(h);
  
  kernel_type Kernel;
  neighsearch_type Nsearch;

  ///
  /// the size of the neighbour list doesn't really
  /// affect the performance of the neighbour search,
  /// so it can be choosen quite large
  ///
  Nsearch.prepare();
  //Nsearch.neighbourList.resize(16384);
  //Nsearch.neighDistList.resize(16384);
  Nsearch.neighbourList.resize(163840);
  Nsearch.neighDistList.resize(163840);
  Logger << "Rankspace prepared";

  ///
  /// 1st SPH sum: density
  /// I need    : pos, h, m
  /// I provide : rho
  ///
  for (size_t k = 0; k < noParts; k++)
    {
      h(k) /= 2.;

      ///
      /// find neighbours
      ///
      const size_t i = k;
      const valueType hi = h(i);
      const valueType srchRad = 2. * hi;
      Nsearch.findNeighbours(i, srchRad);

      const size_t noNeighs = Nsearch.neighbourList[0];
      noneigh(i) = noNeighs;

      static valueType rhoi;
      rhoi = 0.;

      ///
      /// sum over neighbours
      ///
      for (size_t curNeigh = 1; curNeigh <= noNeighs; curNeigh++)
        {
          const valueType r = Nsearch.neighDistList[curNeigh];
          const size_t j = Nsearch.neighbourList[curNeigh];

          //const valueType hij = 0.5 * (hi + h(j));
          const valueType hij = hi;

          rhoi += Kernel.value(r, hij);
        }
      rho(i) = rhoi;
      std::cout << std::scientific << std::setprecision(10)
                << pos(i, X) << " " << pos(i, Y) << " "
                << h(i) << " " << rho(i) << " "
                << std::fixed << std::setw(10) 
                << id(i) << " " << noneigh(i) << "\n";
    }
  Logger << "SPH sum: rho";

  ///
  /// store dump
  ///
  quantsType saveQuants;
  saveQuants.vects += &pos;
  saveQuants.scalars += &rho, &h;
  saveQuants.ints += &id, &noneigh;

  IOManager.saveDump("no_density.h5part", saveQuants);

  MPI::Finalize();
  return EXIT_SUCCESS;
}



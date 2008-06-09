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
typedef sphlatch::valueType valueType;
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

#include "integrator_verlet.h"

#include <boost/progress.hpp>
#include <vector>

// tree stuff
#include "bhtree.h"

void derivate()
{
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());
  
  //matrixRefType pos(PartManager.pos);
  //matrixRefType vel(PartManager.vel);
  //matrixRefType acc(PartManager.acc);
  
  valvectRefType m(PartManager.m);
  
  idvectRefType id(PartManager.id);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  //const size_t myDomain = CommManager.getMyDomain();
  
  //const valueType gravTheta = 0.7;
  const valueType gravTheta = 0.7;
  const valueType gravConst = 6.674e-11;

//sphlatch::BHtree<sphlatch::Monopoles>   Tree(gravTheta,
  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
//sphlatch::BHtree<sphlatch::Octupoles>   Tree(gravTheta,
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

  Tree.detParticleOrder();
  Tree.calcMultipoles();
  
  ///
  /// calculate the accelerations
  ///
  for (size_t i = 0; i < noParts; i++)
  {
    const size_t curIndex = Tree.particleOrder[i];
    ///
    /// don't calculate the acceleration for particle with ID = 1 (star)
    ///
    if ( lrint( id(curIndex) ) != 1 )
    {
      Tree.calcGravity( curIndex );
    }
  }
};

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
  
  size_t& step(PartManager.step);
  
  const size_t myDomain = CommManager.getMyDomain();

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos, &vel;
  CommManager.exchangeQuants.scalars += &eps, &m;
  CommManager.exchangeQuants.ints += &id;

  ///
  /// instantate the MetaIntegrator and
  /// register spatial integration
  ///
  sphlatch::VerletMetaIntegrator Integrator(derivate);
  Integrator.regIntegration(pos, vel, acc);

  IOManager.loadDump("nbody_small.h5part");

  ///
  /// define the quantities to save in a dump
  ///
  quantsType saveQuants;
  saveQuants.vects   += &pos, &vel;
  saveQuants.scalars += &eps, &m;
  saveQuants.ints    += &id;

  ///
  /// exchange particles
  ///
  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());

  ///
  /// prepare ghost sends and send ghosts
  ///
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);
  CommManager.sendGhosts(m);

  ///
  /// define some simulation parameters
  /// 
  const valueType dt = 138.e6 / 50.;
  //const size_t maxSteps = 500000;
  //const size_t saveSteps =  5000;
  const size_t maxStep = 100;
  const size_t saveSteps =  10;

  ///
  /// bootstrap the integrator
  ///
  Integrator.bootstrap(dt);

  ///
  /// the integration loop
  ///
  while ( step < maxStep )
  {
    ///
    /// exchange particles and ghosts
    ///
    CostZone.createDomainPartsIndex();
    CommManager.exchange(CostZone.domainPartsIndex,
                         CostZone.getNoGhosts());
    
    CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  
    CommManager.sendGhosts(pos);
    CommManager.sendGhosts(id);
    CommManager.sendGhosts(m);
   
    ///
    /// integrate
    ///
    Integrator.integrate(dt);
    step++;

    ///
    /// save a dump
    ///
    if ( ( step % saveSteps ) == 0 )
    {
      std::cerr << "save dump!\n";
      IOManager.saveDump("saveDump.h5part", saveQuants);
    }
  }


  MPI::Finalize();
  return EXIT_SUCCESS;
}



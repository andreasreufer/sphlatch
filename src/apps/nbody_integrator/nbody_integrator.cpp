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
typedef sphlatch::valueType valueType;
typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::partsIndexVectType partsIndexVectType;

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
  matrixRefType acc(PartManager.acc);
  
  valvectRefType m(PartManager.m);
  
  idvectRefType id(PartManager.id);

  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = noParts + PartManager.getNoGhostParts();
  const size_t myDomain = CommManager.getMyDomain();
  
  //const valueType gravTheta = 0.7;
  const valueType gravTheta = 0.7;
  //const valueType gravTheta = 0.49;
  const valueType gravConst = 6.674e-11;

  sphlatch::BHtree<sphlatch::Quadrupoles> Tree(gravTheta,
                                               gravConst,
                                               CostZone.getDepth(),
                                               CostZone.getCenter(),
                                               CostZone.getSidelength()
                                               );


  for (size_t i = 0; i < noTotParts; i++)
  {
    Tree.insertParticle(i);
  }

  Tree.detParticleOrder();
  Tree.calcMultipoles();
  
  for (size_t i = 0; i < noParts; i++)
  {
    const size_t curIndex = Tree.particleOrder[i];
    if ( lrint( id(curIndex) ) != 1 )
    {
      Tree.calcGravity( curIndex );
    }
    else
    {
      //std::cerr << "star on domain " << myDomain << "\n";
    }
    if ( lrint( id(curIndex) ) == 3 )
    {
      //std::cout << acc(curIndex, 0) << " " << acc(curIndex, 1) << " " << acc(curIndex, 2) << "\n";
    }
  }
    
  /*
  const valueType dist = sqrt( pos(0,0)* pos(0,0) +  pos(0,1)*pos(0,1));

  std::cerr << " part dist is " << dist << "\n";
  
  const valueType distPow3 = dist*dist*dist;

  std::cerr << " dev: " << ( acc(0, 0) + ( pos(0,0) / distPow3 ) )
            <<    "   " << ( acc(0, 1) + ( pos(0,1) / distPow3 ) ) << "\n";

  std::cerr << " acc: " << acc(0, 0)
            <<    "   " << acc(0, 1) << "\n";

  std::cerr << "\n\n";*/
  /*
  acc(1, 0) = 0.;
  acc(1, 1) = 0.;
  acc(1, 2) = 0.;
  */

    //std::cout << "acc " << acc(0, sphlatch::X) << "\t" << acc(0, sphlatch::Y) << "\n";

  /*for (size_t i = 0; i < noParts; i++)
  {
    const valueType dist = sqrt( pos(i, sphlatch::X)*pos(i, sphlatch::X) +
                                 pos(i, sphlatch::Y)*pos(i, sphlatch::Y) );
    if ( dist > 1.e-9 )
    {
    const valueType distPow3 = dist*dist*dist;
    
    acc(i, sphlatch::X) = - pos(i, sphlatch::X) / distPow3;
    acc(i, sphlatch::Y) = - pos(i, sphlatch::Y) / distPow3;
    acc(i, sphlatch::Z) = 0;
    }
  }*/

};

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  using namespace boost::assign;

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

  //IOManager.loadDump("nbody_start.h5part");
  IOManager.loadDump("nbody_small.h5part");
  //size_t noParts = 200;
  /*size_t noParts = 5;

  //PartManager.setNoParts(1000000,1000000);
  PartManager.setNoParts(noParts,10);
  //PartManager.setNoParts(noParts);
  PartManager.resizeAll();

  //const valueType k = 1.88;
  const valueType k = 0.;
  for (size_t i = 0; i < noParts; i++)
  {
    const valueType phi = 2*M_PI*static_cast<valueType>(random())/RAND_MAX;

    m(i) = 1.e-1;
    eps(i) = 0.;
    id(i) = i;

    pos(i, sphlatch::X) = sin( phi );
    pos(i, sphlatch::Y) = cos( phi );
    pos(i, sphlatch::Z) = 0;

    const valueType theta = phi
      + k*( (static_cast<valueType>(random())/RAND_MAX) - 0.5);
    
    //vel(i, X) = 1.40*cos( theta );
    //vel(i, Y) = -1.35*sin( theta );
    vel(i, sphlatch::X) = cos( theta );
    vel(i, sphlatch::Y) = -sin( theta );
    vel(i, sphlatch::Z) = 0.;
    //std::cerr << phi << " " << theta << "\n";
  }
 
  pos(1, 0) = 0.;
  pos(1, 1) = 0.;
  pos(1, 2) = 0.;
  vel(1, 0) = 0.;
  vel(1, 1) = 0.;
  vel(1, 2) = 0.;
  m(1) = 1.;
  */

  const size_t myDomain = CommManager.getMyDomain();

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
  
  //valueType dt = 2*M_PI / 50.;
  valueType dt = 138.e6 / 50.;

  Integrator.bootstrap(dt);


  while ( step < 10000 )
  {
    CostZone.createDomainPartsIndex();
    CommManager.exchange(CostZone.domainPartsIndex,
                         CostZone.getNoGhosts());
    
    CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  
    CommManager.sendGhosts(pos);
    CommManager.sendGhosts(id);
    CommManager.sendGhosts(m);
    

    for (size_t i = 0; i < PartManager.getNoLocalParts(); i++)
    {
      if ( lrint( id(i) ) == 3 )
      {
        std::cout << pos(i, 0) << " " << pos(i, 1) << " " << pos(i, 2) << " #pos\n";
        std::cout << vel(i, 0) << " " << vel(i, 1) << " " << vel(i, 2) << " #vel\n";
        std::cout << "#" << step << "\n";
        //std::cout << vel(i, 0) << " " << vel(i, 1) << " " << vel(i, 2) << "\n";

        const valueType dist = sqrt( pos(i, 0)*pos(i, 0) +
                                     pos(i, 1)*pos(i, 1) +
                                     pos(i, 2)*pos(i, 2) );
        const valueType vel1 = sqrt( vel(i, 0)*vel(i, 0) +
                                     vel(i, 1)*vel(i, 1) +
                                     vel(i, 2)*vel(i, 2) );

        //std::cout << dist << " " << vel1 << "     " <<  pos(i, 0) << "\n";//   (" << myDomain << ")\n";
      }
    }
    
    Integrator.integrate(dt);
    
    step++;
  }


  MPI::Finalize();
  return EXIT_SUCCESS;
}



// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

// enable parallel version
//#define SPHLATCH_PARALLEL

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
//typedef sphlatch::fType fType;
//typedef sphlatch::partsIndexVectType partsIndexVectType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

//#include "communication_manager.h"
//typedef sphlatch::CommunicationManager com_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include "integrator_verlet.h"
#include "integrator_predcorr.h"

#include <boost/progress.hpp>
#include <vector>

// tree stuff
#include "bhtree.h"

//#define ORBIT
#define ONEDIMTHROW

void derivate()
{
  using namespace sphlatch;
  
  part_type& PartManager(part_type::instance());
  matrixRefType pos(PartManager.pos);
  //matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  const size_t noParts = PartManager.getNoLocalParts();

  for (size_t i = 0; i < noParts; i++)
  {
#ifdef ORBIT
    const fType dist = sqrt( pos(i, X)*pos(i, X) +
                                 pos(i, Y)*pos(i, Y) );
    
    const fType distPow3 = dist*dist*dist;
    
    acc(i, X) = - pos(i, X) / distPow3;
    acc(i, Y) = - pos(i, Y) / distPow3;
    acc(i, Z) = 0;
#endif
#ifdef ONEDIMTHROW
    //const fType dist = sqrt( pos(i, X)*pos(i, X) );
    //const fType distPow3 = dist*dist*dist;

    //acc(i, X) = - pos(i, X) / distPow3;
    acc(i, X) = (pos(i, X) > 0.) ? -1. : 1.;
#endif
  }

};

int main(int argc, char* argv[])
{

  using namespace boost::assign;
  using namespace sphlatch;

  part_type& PartManager(part_type::instance());

  PartManager.useGravity();
  PartManager.setNoParts(100);
  PartManager.resizeAll();
  

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType acc(PartManager.acc);

  //valvectRefType m(PartManager.m);

  idvectRefType id(PartManager.id);

  //VerletMetaIntegrator myIntegrator(derivate);
  PredCorrMetaIntegrator myIntegrator(derivate);
  myIntegrator.regIntegration(pos, vel, acc);

  //size_t noParts = 200;
  size_t noParts = 1;

  //PartManager.setNoParts(1000000,1000000);
  PartManager.setNoParts(noParts,10);
  PartManager.resizeAll();

#ifdef ORBIT
  //const fType k = 1.88;
  const fType k = 0.;
  for (size_t i = 0; i < noParts; i++)
  {
    id(i) = i;

    const fType phi = 2*M_PI*static_cast<fType>(random())/RAND_MAX;

    pos(i, X) = sin( phi );
    pos(i, Y) = cos( phi );
    pos(i, Z) = 0;

    const fType theta = phi
      + k*( (static_cast<fType>(random())/RAND_MAX) - 0.5);
    
    //vel(i, X) = 1.40*cos( theta );
    //vel(i, Y) = -1.35*sin( theta );
    vel(i, X) = cos( theta );
    vel(i, Y) = -sin( theta );
    vel(i, Z) = 0.;
    std::cerr << phi << " " << theta << "\n";
  }
#endif
#ifdef ONEDIMTHROW
  for (size_t i = 0; i < noParts; i++)
  {
    id(i) = i;

    pos(i, X) = 1.;
    pos(i, Y) = 0.;
    pos(i, Z) = 0.;

    vel(i, X) = 0.;
    vel(i, Y) = 0.;
    vel(i, Z) = 0.;
  }
#endif
  
  const fType dt = 2*M_PI / 50.;
  
  myIntegrator.bootstrap(dt);

  size_t step = 0;
  while ( step < 1000 )
  {
    myIntegrator.integrate(dt);
    step++;
    
    for (size_t i = 0; i < noParts; i++)
    {
      std::cout << pos(i, X) << "\t" << pos(i, Y) << "\n";
    }
  }

  std::cerr << pos << "\n";
  std::cerr << vel << "\n";
  std::cerr << acc << "\n";

  return EXIT_SUCCESS;
}



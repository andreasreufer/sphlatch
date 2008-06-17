// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#include <iostream>
#include <iomanip>

#include "typedefs.h"
//typedef sphlatch::valueType valueType;
//typedef sphlatch::partsIndexVectType partsIndexVectType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "kernel_cubicspline3d.h"

int main(int argc, char* argv[])
{
  using namespace boost::assign;
  using namespace sphlatch;

  part_type& PartManager(part_type::instance());

  PartManager.useBasicSPH();
  PartManager.setNoParts(100);
  PartManager.resizeAll();

  CubicSpline3D myKernel;
  
  matrixRefType pos(PartManager.pos);
  valvectRefType h(PartManager.m);

  valvectType rvec(3);
  rvec(X) = 0.;
  rvec(Y) = 0.;
  rvec(Z) = 0.;

  for (size_t j = 0; j < 10; j++)
  {
  //const size_t steps = 10000;
  const size_t steps = 10000000;
  for (size_t i = 0; i < steps; i++)
  {
    const valueType x = -2.25 + ( 4.5*i / steps );
    const valueType r = fabs(x);
    const valueType h = 1.;
    
    rvec(Y) = x;
    
    //std::cout << x << "\t" << myKernel.value( r, h ) << "\t" << myKernel.derive( r, h, rvec )(Y) << "\n";
    myKernel.value( r, h );
    myKernel.derive( r, h, rvec );
  }
  }

  return EXIT_SUCCESS;
}



// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "kernel_cubicspline3d.h"

enum vectIndices { X, Y, Z };

int main(int argc, char* argv[])
{
  using namespace sphlatch;

  CubicSpline3D myKernel;
  
  valvectType rvec(3);
  rvec(X) = 0.;
  rvec(Y) = 0.;
  rvec(Z) = 0.;

  const size_t steps = 1000;
  const valueType h = 1.e-9;
  
  valueType kernelVal = 0;
  valvectType derivVal(3);
  derivVal(X) = 0.;
  derivVal(Y) = 0.;
  derivVal(Z) = 0.;
  
  /*for (size_t i = 0; i < steps; i++)
  {
    const valueType x = 4.5*( static_cast<valueType>( rand() ) / RAND_MAX ) - 2.25; 
    const valueType y = 4.5*( static_cast<valueType>( rand() ) / RAND_MAX ) - 2.25; 
    const valueType z = 4.5*( static_cast<valueType>( rand() ) / RAND_MAX ) - 2.25; 
    
    const valueType r = sqrt( x*x + y*y + z*z );
  
    rvec(X) = x;
    rvec(Y) = y;
    rvec(Z) = z;

    kernelVal += myKernel.value( r, h )*4.5*4.5*4.5;
    derivVal  += myKernel.derive( r, h, rvec )*4.5*4.5*4.5;
  
  }
  kernelVal /= static_cast<valueType>( steps );
  derivVal  /= static_cast<valueType>( steps );
  */
  
  const valueType dx = 4.5*h / static_cast<valueType>(steps);
  const valueType vol = 4.5*4.5*4.5*h*h*h;
  
  for (size_t i = 0; i < steps; i++)
  {
    const valueType x = - 2.25*h + ( 0.5 + i )*dx;
    for (size_t j = 0; j < steps; j++)
    {
      const valueType y = - 2.25*h + ( 0.5 + j )*dx;
      for (size_t k = 0; k < steps; k++)
      {
        const valueType z = - 2.25*h + ( 0.5 + k )*dx;

        const valueType r = sqrt( x*x + y*y + z*z );
  
        rvec(X) = x;
        rvec(Y) = y;
        rvec(Z) = z;

        kernelVal += myKernel.value( r, h )*vol;
        derivVal  += myKernel.derive( r, h, rvec )*vol;
      }
    }
    std::cout << i << " " << std::flush;
  }
  std::cout << "\n";

  kernelVal /= static_cast<valueType>( steps*steps*steps );
  derivVal  /= static_cast<valueType>( steps*steps*steps );

  std::cout << "kernel       integral: " << kernelVal << "\n";
  std::cout << "kernel deriv integral: " << derivVal << "\n";

  return EXIT_SUCCESS;
}



// some defs

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

#include <iostream>
#include <iomanip>

#include "typedefs.h"

#include "kernel_cubicspline.h"

enum vectIndices { X, Y, Z };

int main(int argc, char* argv[])
{
  using namespace sphlatch;

  CubicSpline3D myKernel;
//CubicSpline2D myKernel;

  const size_t steps = 200;
  const fType h = 17.;

  fType kernelVal = 0;

  fType derivValX = 0., derivValY = 0., derivValZ = 0.;

  const fType dx = 4.5 * h / static_cast<fType>(steps);
  const fType vol = dx * dx * dx;
  //const fType vol = dx * dx;

  for (size_t i = 0; i < steps; i++)
    {
      const fType x = -2.25 * h + (0.5 + i) * dx;
      for (size_t j = 0; j < steps; j++)
        {
          const fType y = -2.25 * h + (0.5 + j) * dx;
          for (size_t k = 0; k < steps; k++)
            {
              const fType z = -2.25 * h + (0.5 + k) * dx;

              const fType r = sqrt(x * x + y * y + z * z);

              kernelVal += myKernel.value(r, h) * vol;

              myKernel.derive(r, h, x, y, z);
              derivValX += myKernel.derivX * vol;
              derivValY += myKernel.derivY * vol;
              derivValZ += myKernel.derivZ * vol;
            }

          /*const fType r = sqrt(x * x + y * y);

          kernelVal += (myKernel.value(r, h) * vol);

          myKernel.derive(r, h, x, y);
          derivValX += myKernel.derivX * vol;
          derivValY += myKernel.derivY * vol;*/
        }
      std::cout << i << " " << std::flush;
    }
  std::cout << "\n";

  std::cout << "kernel       integral: " << kernelVal << "\n";
  std::cout << "kernel deriv integral: " << derivValX << " " << derivValY << " " << derivValZ << "\n";


  const fType r = 3.0;
  const fType x = 2.0000;
  const fType y = 1.4142;
  const fType z = 1.7321;
  //myKernel.derive(r, h, x, y, z);
  //std::cout << "kernel deriv: [" << myKernel.derivX << "," << myKernel.derivY << "," << myKernel.derivZ << "]\n";

  return EXIT_SUCCESS;
}



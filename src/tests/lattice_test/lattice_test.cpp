// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#include <iostream>
#include <iomanip>

#include "typedefs.h"
typedef sphlatch::valueType valueType;
typedef sphlatch::valvectType valvectType;

#include "lattice_hcp.h"
typedef sphlatch::LatticeHCP lattice_type;

//using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  lattice_type Lattice(1., 30., 30., 30.);

  size_t partCount = 0;

  while (!Lattice.isLast)
    {
      if (Lattice.rCur < 20. && Lattice.zCur < 0.1)
        {
          partCount++;
          std::cout << Lattice.xCur << "\t"
                    << Lattice.yCur << "\t"
                    << Lattice.zCur << "\n";
        }
      Lattice.next();
    }
  std::cerr << partCount << " particles!\n";

  return EXIT_SUCCESS;
}



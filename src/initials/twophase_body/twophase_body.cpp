// some defs

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/std/vector.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::valueType valueType;
typedef sphlatch::valvectType valvectType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "lagrange_sphere1D_solver.h"
typedef sphlatch::LagrangeSphere1DSolver lg1D_solver_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options()
  ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), "output file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type&        IOManager(io_type::instance());
  part_type&      PartManager(part_type::instance());

  const valueType mTot = 5.3014e27;
  const valueType mCore = 0.30 * mTot;
  const valueType mMantle = 0.70 * mTot;

  ///
  /// guess for the density
  ///
  const valueType rhoCore = 6.;
  const valueType rhoMantle = 3.;

  const valueType rCore = pow(0.75 * mCore / (M_PI * rhoCore), 1. / 3.);
  const valueType rMantle = pow(0.75 * mMantle / (M_PI * rhoMantle)
                                + pow(rCore, 3.), 1. / 3.);
  const size_t noCells = 200;
  const size_t noEdges = noCells + 1;

  const size_t noCoreCells = lrint(noCells * (rCore / rMantle));
  const size_t noMantleCells = noCells - noCoreCells;

  lg1D_solver_type Solver(noCells);

  const valueType drCore = rCore / noCoreCells;
  const valueType drMantle = (rMantle - rCore) / noMantleCells;

  ///
  /// set the shell edges
  ///
  Solver.r(0) = 0.;
  Solver.v(0) = 0.;
  for (size_t i = 1; i < noEdges; i++)
    {
      if (i <= noCoreCells)
        Solver.r(i) = Solver.r(i - 1) + drCore;
      else
        Solver.r(i) = Solver.r(i - 1) + drMantle;

      Solver.v(i) = 0.;
    }

  ///
  /// set the shell centered values
  ///
  for (size_t i = 0; i < noCells; i++)
    {
      const valueType rI = Solver.r(i);
      const valueType rO = Solver.r(i + 1);
      const valueType rC = 0.5 * (rI + rO);

      const valueType vol = (4. * M_PI / 3.) * (rO * rO * rO - rI * rI * rI);

      if (rC < rCore)
        {
          Solver.m(i) = vol * rhoCore;
          Solver.mat(i) = 5;
        }
      else
        {
          Solver.m(i) = vol * rhoMantle;
          Solver.mat(i) = 4;
        }

      Solver.u(i) = 1.e10;
    }

  ///
  /// set some constants
  ///
  Solver.uMin = 1.e4;
  Solver.friction = 1.e-2;

  ///
  /// integrate for a certain physical time
  ///
  Solver.integrateTo(5.e3);

  valvectType rCen;
  rCen.resize(noCells);

  for (size_t i = 0; i < noCells; i++)
    {
      rCen(i) = 0.5 * (Solver.r(i) + Solver.r(i + 1));
    }

  std::string dumpfilename = "dump.hdf5";
  IOManager.savePrimitive(Solver.r, "r", dumpfilename);
  IOManager.savePrimitive(rCen, "rCen", dumpfilename);
  IOManager.savePrimitive(Solver.v, "v", dumpfilename);
  IOManager.savePrimitive(Solver.m, "m", dumpfilename);
  IOManager.savePrimitive(Solver.u, "u", dumpfilename);
  IOManager.savePrimitive(Solver.p, "p", dumpfilename);
  IOManager.savePrimitive(Solver.rho, "rho", dumpfilename);


  const valueType sq32 = sqrt(3.) / 2.;
  const valueType sq34 = sqrt(3.) / 4.;
  const valueType rMax = 100.;

  valueType xCur = 0, yCur = 0, zCur = 0;
  int ix = 0, iy = 0, iz = 0;

  ///
  /// go to lattice point with the smallest coordinates
  ///
  while (xCur > -rMax)
    {
      xCur -= 1.;
      ix -= 1;
    }

  while (yCur > -rMax)
    {
      yCur -= sq32;
      iy -= 1;
    }
  if (iy % 2 == 1)
    xCur -= 0.5;

  while (zCur > -rMax)
    {
      zCur -= sq34;
      iz -= 1;
    }
  if (iz % 2 == 1)
    {
      xCur -= 0.5;
      yCur -= sq34;
    }

  const int noIx = 2 * (-ix);
  const int noIy = 2 * (-iy);
  const int noIz = 2 * (-iz);

  const int ixMax = -ix;
  const int iyMax = -iy;
  const int izMax = -iz;

  size_t countParts = 0;
  while (iz <= izMax)
    {
      const int prevIy = iy;
      const valueType yPrev = yCur;

      while (iy <= iyMax)
        {
          const int prevIx = ix;
          const valueType xPrev = xCur;
          while (ix <= ixMax)
            {
              const valueType rCur = sqrt(xCur * xCur +
                                          yCur * yCur +
                                          zCur * zCur);
              if (rCur < rMax)
                {
                  std::cout << xCur << "\t"
                            << yCur << "\t"
                            << zCur << "\n";
                  countParts++;
                }
              xCur += 1.;
              ix++;
            }
          xCur = xPrev;
          ix = prevIx;

          iy++;
          if (iy <= iyMax)
            {
              if ((iy % 2) == 0)
                {
                  xCur += 0.5;
                }
              else
                {
                  xCur -= 0.5;
                }
              yCur += sq32;
            }
        }
      yCur = yPrev;
      iy = prevIy;

      iz++;
      zCur += sq34;
      if (iy <= iyMax)
        {
          if ((iz % 2) == 0)
            {
              xCur += 0.5;
              yCur += sq34;
            }
          else
            {
              xCur -= 0.5;
              yCur -= sq34;
            }
        }
    }

  std::cerr << countParts << "\n";


  /*sphlatch::quantsType saveQuants;
     PartManager.step = 0;

     std::string outputFilename = VMap["output-file"].as<std::string>();
     std::cout << " -> " << outputFilename << "\n";
     IOManager.saveDump( outputFilename, saveQuants );*/

  return EXIT_SUCCESS;
}



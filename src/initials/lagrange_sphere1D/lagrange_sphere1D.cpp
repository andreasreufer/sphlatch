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

  const size_t noCells = 400;
  const valueType dr = 6.e8 / noCells;

  lg1D_solver_type Solver(noCells);

  Solver.uMin = 1.e4;
  //Solver.gravConst = 0.;
  Solver.friction = 1.e-2;

  for (size_t i = 0; i < noCells; i++)
  {
    Solver.r(i) = dr*(i);
    Solver.v(i) = 0.;
   
    if ( i < 200 )
    {
      Solver.rho(i) = 6.;
      Solver.u(i)   = 1.e10;
      Solver.mat(i) = 5;
    }
    else
    {
      Solver.rho(i) = 3.;
      Solver.u(i)   = 1.e10;
      Solver.mat(i) = 4;
    }

  }
  Solver.r(0) = 0.;
  Solver.r(noCells) = (noCells)*dr;
  Solver.v(noCells) = 0.;
  
  /*for (size_t i = 400; i < 600; i++)
  {
    Solver.u(i) = 7.0e10;
  }*/
  
  Solver.densityToMass();

  std::string dumpfilename;
  dumpfilename = "dump01.hdf5";
  IOManager.savePrimitive( Solver.r, "r", dumpfilename);
  IOManager.savePrimitive( Solver.v, "v", dumpfilename);
  IOManager.savePrimitive( Solver.m, "m", dumpfilename);
  IOManager.savePrimitive( Solver.u, "u", dumpfilename);
  IOManager.savePrimitive( Solver.p, "p", dumpfilename);
  IOManager.savePrimitive( Solver.rho, "rho", dumpfilename);
  
  Solver.integrateTo(5.e3);
  
  dumpfilename = "dump02.hdf5";
  IOManager.savePrimitive( Solver.r, "r", dumpfilename);
  IOManager.savePrimitive( Solver.v, "v", dumpfilename);
  IOManager.savePrimitive( Solver.m, "m", dumpfilename);
  IOManager.savePrimitive( Solver.u, "u", dumpfilename);
  IOManager.savePrimitive( Solver.p, "p", dumpfilename);
  IOManager.savePrimitive( Solver.rho, "rho", dumpfilename);

  sphlatch::quantsType saveQuants;
  PartManager.step = 0;
  
  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cout << " -> " << outputFilename << "\n";
  IOManager.saveDump( outputFilename, saveQuants );

  return EXIT_SUCCESS;
}



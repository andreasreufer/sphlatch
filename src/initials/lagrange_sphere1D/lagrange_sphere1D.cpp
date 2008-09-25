// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
//#define SPHLATCH_PARALLEL

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
typedef sphlatch::identType identType;
typedef sphlatch::valvectType valvectType;

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "lagrange_sphere1D_solver.h"
typedef sphlatch::LagrangeSphere1DSolver lg1D_solver_type;

#ifdef SPHLATCH_ANEOS
#include "eos_aneos.h"
typedef sphlatch::ANEOS eosType;
#endif

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  //MPI::Init(argc, argv);

  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), " output file")
  ("mtot,M", po::value<valueType>(),          " total  mass")
  ("core-frac,f", po::value<valueType>(),     " core   mass fraction")
  ("mat-core", po::value<identType>(),        " core   material ")
  ("mat-mantle", po::value<identType>(),      " mantle material")
  ("rho-core", po::value<valueType>(),        " core   initial density")
  ("rho-mantle", po::value<valueType>(),      " mantle initial density ")
#ifdef SPHLATCH_ANEOS
  ("T",          po::value<valueType>(),      " initial temperature in eV")
#else
  ("u",          po::value<valueType>(),      " initial specific energy")
#endif
  ("umin",       po::value<valueType>(),      " minimal specific energy")
  ("fric",       po::value<valueType>(),      " inverse friction timescale");

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

  valueType mTot = 5.3014e27;
  if (VMap.count("mtot"))
    mTot = VMap["mtot"].as<valueType>();
  
  valueType mCoreFrac = 0.30;
  if (VMap.count("core-frac"))
    mCoreFrac = VMap["core-frac"].as<valueType>();
  
  const valueType mCore  = mCoreFrac * mTot;

  identType matCore = 5;
  if (VMap.count("mat-core"))
    matCore = VMap["mat-core"].as<identType>();
  
  identType matMantle = 4;
  if (VMap.count("mat-mantle"))
    matMantle = VMap["mat-mantle"].as<identType>();

  valueType rhoCore = 6.;
  if (VMap.count("rho-core"))
    rhoCore = VMap["rho-core"].as<valueType>();
  
  valueType rhoMantle = 3.;
  if (VMap.count("rho-mantle"))
    rhoMantle = VMap["rho-mantle"].as<valueType>();
 
#ifdef SPHLATCH_ANEOS
  valueType TInit = 0.1;
  if (VMap.count("T"))
    TInit = VMap["T"].as<valueType>();
#else 
  valueType uInit = 5.e10;
  if (VMap.count("u"))
    uInit = VMap["u"].as<valueType>();
#endif

  valueType uMin = 1.e4;
  if (VMap.count("umin"))
    uMin = VMap["umin"].as<valueType>();

  valueType fric = 1.e-2;
  if (VMap.count("fric"))
    fric = VMap["fric"].as<valueType>();

  ///
  /// guess for the density
  ///
#ifdef SPHLATCH_ANEOS
  const size_t noCells = 300;
#else
  const size_t noCells = 1000;
#endif
  const size_t noEdges = noCells + 1;
  const size_t noCoreCells = lrint(noCells * (mCore / mTot));
  const valueType dm = mTot / static_cast<valueType>(noCells);

  std::cerr << " total  mass:          " << mTot << "\n"
            << " core   mass:          " << mCore << "   ("
            << noCoreCells << " cells)\n"
            << " core   material:      " << matCore << "\n"
            << " mantle mass:          " << mTot - mCore << "   ("
            << noCells - noCoreCells << " cells)\n"
            << " mantle material:      " << matMantle << "\n"
#ifdef SPHLATCH_ANEOS
            << " initial T:            " << TInit     << "\n"
#else
            << " initial u:            " << uInit     << "\n"
#endif
            << " minimal u:            " << uMin      << "\n"
            << " inv. fric. timescale: " << fric << "\n\n";

  ///
  /// instantate solver
  ///
  lg1D_solver_type Solver(noCells);

  ///
  /// set the shell edges
  ///
  Solver.r(0) = 0.;
  Solver.v(0) = 0.;
  for (size_t i = 1; i < noEdges; i++)
    {
      valueType rhoCur = 0.;
      if (i <= noCoreCells)
        rhoCur = rhoCore;
      else
        rhoCur = rhoMantle;
      Solver.r(i) = pow( pow(Solver.r(i-1),3.) +
                         (3./(4*M_PI))*(dm / rhoCur), 1./3.);
      Solver.v(i) = 0.;
    }
  
  ///
  /// set the shell centered values
  ///
#ifdef SPHLATCH_ANEOS
  eosType EOS(eosType::instance());
#endif
for (size_t i = 0; i < noCells-1; i++)
    {
      Solver.m(i) = dm;

      if (i < noCoreCells)
        {
          Solver.mat(i) = matCore;
          Solver.rho(i) = rhoCore;
        }
      else
        {
          Solver.mat(i) = matMantle;
          Solver.rho(i) = rhoMantle;
        }
#ifdef SPHLATCH_ANEOS
      static valueType curP, curCs;
      EOS.getSpecEnergy( Solver.rho(i), TInit, Solver.mat(i),
                         curP, curCs, Solver.u(i) );
#else
      Solver.u(i) = uInit;
#endif
    }
  Solver.m(noCells-1) = 0.;
  Solver.u(noCells-1) = 0.;
  Solver.mat(noCells-1) = matMantle;

  ///
  /// set some constants
  ///
  Solver.uMin = uMin;
  Solver.friction = fric;

  ///
  /// integrate for a certain physical time
  ///
  std::cerr << " start 1D Lagrange solver\n";
  Solver.integrateTo(5.e3); /// replace by option
  //Solver.integrateTo(1.e-1); /// replace by option
  std::cerr << " ... finished\n";
  
  ///
  /// the last cell contains vacuum, which will assign
  /// far outside zero mass. circumvent this, by setting
  /// the same density and temperature like on the second
  /// to last cell
  ///
  Solver.rho(noCells-1) = Solver.rho(noCells-2);
  Solver.u(noCells-1) = Solver.u(noCells-2);

  ///
  /// calculate the position of the cell centers
  ///
  valvectType rCen;
  rCen.resize(noCells);

  for (size_t i = 0; i < noCells; i++)
    {
      rCen(i) = 0.5 * (Solver.r(i) + Solver.r(i + 1));
    }

  std::string dumpfilename = VMap["output-file"].as<std::string>();
  IOManager.savePrimitive(rCen, "r", dumpfilename);
  IOManager.savePrimitive(Solver.v, "v", dumpfilename);
  IOManager.savePrimitive(Solver.m, "m", dumpfilename);
  IOManager.savePrimitive(Solver.u, "u", dumpfilename);
  IOManager.savePrimitive(Solver.p, "p", dumpfilename);
  IOManager.savePrimitive(Solver.mat, "mat", dumpfilename);
  IOManager.savePrimitive(Solver.rho, "rho", dumpfilename);
#ifdef SPHLATCH_ANEOS
  IOManager.savePrimitive(Solver.T, "T", dumpfilename);
  IOManager.savePrimitive(Solver.phase, "phase", dumpfilename);
#endif
  IOManager.saveAttribute(Solver.gravConst, "gravconst", dumpfilename);
  IOManager.saveAttribute(Solver.uMin, "umin", dumpfilename);

  //MPI::Finalize();
  return EXIT_SUCCESS;
}


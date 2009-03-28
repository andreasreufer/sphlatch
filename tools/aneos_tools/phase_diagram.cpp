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

#include <boost/progress.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::fType fType;
typedef sphlatch::idType idType;
typedef sphlatch::valvectType valvectType;

#include "eos_aneos.h"
typedef sphlatch::ANEOS eos_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  /*if (VMap.count("help") ||
      not VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }*/

  eos_type&      EOS(eos_type::instance());

  // in:  rho, u mat
  // out: P, cs, T, phase
  //

  const fType evToK = 11605.333;

  const fType rhoMin = 0.5;
  const fType rhoMax = 1.5;
  const size_t rhoSteps = 400;
  
  const fType Tmin = 1.e2 / evToK;
  const fType Tmax = 1.e3 / evToK;
  const size_t Tsteps = 250;

  const idType mat = 2;
  
  const fType rhoLdelta = ( log(rhoMax) - log(rhoMin) ) / rhoSteps;
  const fType TLdelta = ( log(Tmax) - log(Tmin) ) / Tsteps;

  for (size_t i = 0; i < rhoSteps; i++)
  {
    for (size_t j = 0; j < Tsteps; j++)
    {
      static fType curP, curCS, curU;
      static idType curPh;

      const fType curRho = rhoMin*exp(rhoLdelta*i);
      const fType curT   =   Tmin*exp(  TLdelta*j);
  
      EOS.getSpecEnergy(curRho, curT, mat, curP, curCS, curU, curPh);

      std::cerr << std::setw(14) << std::setprecision(6) << std::scientific 
                << curRho << "\t"
                << curT * evToK << "\t"
                << curP << "\t"
                << curCS << "\t"
                << curU << "\t"
                << curPh << "\n";
    }
  }


  /*boost::progress_display partProg(noParts);
  for (size_t i = 0; i < noParts; i++)
    {
      EOS(i, p(i), cs(i));
      ++partProg;
    }

  sphlatch::quantsType saveQuants;

  saveQuants.ints += &mat, &phase;
  saveQuants.scalars += &p, &T, &cs;

  IOManager.saveDump(inputFilename, saveQuants);

  std::cerr << " -> " << inputFilename << "\n";*/

  return EXIT_SUCCESS;
}


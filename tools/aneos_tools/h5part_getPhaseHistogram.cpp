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
typedef sphlatch::valvectType valvectType;
typedef sphlatch::zerovalvectType zerovalvectType;

typedef sphlatch::identType identType;

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::quantsType quantsType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file")
  ("mat-number,n", po::value<identType>(),   "material number");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help") ||
      not VMap.count("input-file") ||
      not VMap.count("output-file") ||
      not VMap.count("mat-number") )
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }


  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  std::string inputFilename = VMap["input-file"].as<std::string>();
  std::string outputFilename = VMap["output-file"].as<std::string>();
  const identType matNumber = VMap["mat-number"].as<identType>();

  PartManager.useBasicSPH();
  PartManager.useMaterials();
  PartManager.usePhase();

  IOManager.loadDump(inputFilename);

  valvectRefType m(PartManager.m);
  
  idvectRefType mat(PartManager.mat);
  idvectRefType phase(PartManager.phase);

  const size_t noParts = PartManager.getNoLocalParts();

  valvectType phaseHistogram = zerovalvectType(6);
  for (size_t i = 0; i < noParts; i++)
    {
      if ( mat(i) == matNumber )
        phaseHistogram( phase(i) - 1 ) += m(i);
    }

  std::fstream fout;
  
  fout.open( outputFilename.c_str(), std::ios::app | std::ios::out);
  fout << std::setw(14) << std::setprecision(6) << std::scientific << PartManager.attributes["TIME"];
  for (size_t i = 0; i < 6; i++)
  {
    fout << std::setw(14) << std::setprecision(6) << std::scientific << phaseHistogram(i);
  }
  fout << "\n";
  fout.close();

  return EXIT_SUCCESS;
}


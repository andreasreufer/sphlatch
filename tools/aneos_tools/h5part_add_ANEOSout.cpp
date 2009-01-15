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

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::quantsType quantsType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

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

  if (VMap.count("help") ||
      not VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }


  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());
  eos_type&      EOS(eos_type::instance());

  std::string inputFilename = VMap["input-file"].as<std::string>();

  std::cerr << "add T, p and phase to " << inputFilename << " ... \n";

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
  PartManager.usePhase();
  PartManager.useTemperature();

  IOManager.loadDump(inputFilename);

  valvectRefType p(PartManager.p);
  valvectRefType T(PartManager.T);
  valvectRefType cs(PartManager.cs);

  idvectRefType mat(PartManager.mat);
  idvectRefType phase(PartManager.phase);

  const size_t noParts = PartManager.getNoLocalParts();

  boost::progress_display partProg(noParts);
  for (size_t i = 0; i < noParts; i++)
  {
    EOS(i, p(i), cs(i) );
    ++partProg;
  }

  sphlatch::quantsType saveQuants;
  
  saveQuants.ints += &mat, &phase;
  saveQuants.scalars += &p, &T, &cs;

  IOManager.saveDump(inputFilename, saveQuants);
  
  std::cerr << " -> " << inputFilename << "\n";

  return EXIT_SUCCESS;
}


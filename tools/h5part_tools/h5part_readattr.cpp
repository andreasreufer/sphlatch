#define SPHLATCH_HDF5

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
#include <boost/assign/std/vector.hpp>
namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::fType fType;

#include "io_particle.h"
class particle :
  public sphlatch::IOPart
{
  public:
  ioVarLT getLoadVars()
  {
    ioVarLT vars;
    return vars;
  }
  
  ioVarLT getSaveVars()
  {
    ioVarLT vars;
    return vars;
  }
};

typedef particle partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

partSetT parts;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options() ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input file")
  ("attr-key,k", po::value<std::string>(), "attribute key");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (!VMap.count("input-file") || !VMap.count("attr-key"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  std::string inputFileName = VMap["input-file"].as<std::string>();
  std::string attrKey = VMap["attr-key"].as<std::string>();

  parts.loadHDF5(inputFileName);

  sphlatch::attrMT::iterator attItr = parts.attributes.find(attrKey);

  if ( attItr != parts.attributes.end() )
  {
    std::cout << attrKey << ":\t" 
              << std::scientific << std::setprecision(12) 
              << (*attItr).second << "\n";
    return EXIT_SUCCESS;
  }
  else
  {
      std::cerr << "attribute \"" << attrKey
                << "\" not found in " << inputFileName << "!\n";
      return EXIT_FAILURE;
  }
}



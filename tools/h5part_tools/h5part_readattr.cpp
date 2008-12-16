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

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

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

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  std::string inputFileName = VMap["input-file"].as<std::string>();
  std::string attrKey = VMap["attr-key"].as<std::string>();

  IOManager.loadDump(inputFileName);

  if (PartManager.attrExists(attrKey))
    {
      std::cout << attrKey << ":\t" << PartManager.attributes[attrKey] << "\n";
      return EXIT_SUCCESS;
    }
  else
    {
      std::cerr << "attribute \"" << attrKey
                << "\" not found in " << inputFileName << "!\n";
      return EXIT_FAILURE;
    }
}



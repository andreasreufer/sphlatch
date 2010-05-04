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
typedef sphlatch::quantsType quantsType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;
using namespace sphlatch;

int main(int argc, char* argv[])
{
  typedef sphlatch::fType fType;
  typedef sphlatch::quantsType quantsType;

  po::options_description Options("Global Options");

  Options.add_options() ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input file")
  ("attr-key,k", po::value<std::string>(), "attribute key")
  ("attr-value,v", po::value<fType>(), "attribute value");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (!VMap.count("input-file")
      || !VMap.count("attr-key") || !VMap.count("attr-value"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  std::string inputFileName = VMap["input-file"].as<std::string>();
  std::string attrKey = VMap["attr-key"].as<std::string>();
  fType attrValue = VMap["attr-value"].as<fType>();
  quantsType saveQuants;

  IOManager.loadDump(inputFileName);
  PartManager.attributes[attrKey] = attrValue;
  IOManager.saveDump(inputFileName, saveQuants);
}



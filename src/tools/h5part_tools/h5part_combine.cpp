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
typedef sphlatch::valueType valueType;

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
  ("inputa-file,a", po::value<std::string>(), "input file A")
  ("inputb-file,b", po::value<std::string>(), "inout file B")
  ("output-file,o", po::value<std::string>(), "output file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (!VMap.count("output-file") || VMap.count("help") ||
      !VMap.count("inputa-file") || !VMap.count("inputb-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  std::string inputAFileName = VMap["inputa-file"].as<std::string>();
  std::string inputBFileName = VMap["inputb-file"].as<std::string>();
  std::string outputFileName = VMap["output-file"].as<std::string>();

  PartManager.useBasicSPH();
  PartManager.useMaterials();
  PartManager.useEnergy();

  using namespace sphlatch;
  using namespace boost::assign;

  idvectRefType id(PartManager.id);
  idvectRefType mat(PartManager.mat);

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);

  valvectRefType m(PartManager.m);
  valvectRefType u(PartManager.u);

  quantsType saveQuants;
  saveQuants.ints += &id, &mat;
  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &m, &u;

  IOManager.loadDump(inputAFileName);
  IOManager.loadDump(inputBFileName);
  IOManager.saveDump(outputFileName, saveQuants);


  return EXIT_SUCCESS;
}



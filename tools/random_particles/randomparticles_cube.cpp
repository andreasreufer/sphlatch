// some defs

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

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
  ("no-parts,n",   po::value<size_t>(),       "no. of particles")
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
  
  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  std::string outputFileName = VMap["output-file"].as<std::string>();
  const size_t noParts       = VMap["no-parts"   ].as<size_t>();

  PartManager.useGravity();

  PartManager.setNoParts(noParts, 0);
  PartManager.resizeAll();

  sphlatch::idvectRefType id( PartManager.id );
  sphlatch::matrixRefType pos( PartManager.pos );
  sphlatch::valvectRefType m( PartManager.m );

  for (size_t i = 0; i < noParts; i++)
  {
    id(i) = i;

    m(i) = 1;

    pos(i, X) = static_cast<sphlatch::valueType>(rand()) / RAND_MAX;
    pos(i, Y) = static_cast<sphlatch::valueType>(rand()) / RAND_MAX;
    pos(i, Z) = static_cast<sphlatch::valueType>(rand()) / RAND_MAX;
  }
  

  using namespace boost::assign;
  sphlatch::quantsType saveQuants;
  saveQuants.vects   += &pos;
  saveQuants.scalars += &m;
  saveQuants.ints    += &id;

  PartManager.step = 0;
  IOManager.setSinglePrecOut();
  IOManager.saveDump( outputFileName, saveQuants );

  PartManager.attributes["time"] = 0.0;

 
  return EXIT_SUCCESS;
}



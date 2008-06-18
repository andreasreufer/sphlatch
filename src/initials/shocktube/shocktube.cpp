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
#include <boost/lexical_cast.hpp>
#include <boost/assign/std/vector.hpp>

namespace po = boost::program_options;

#include "typedefs.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;
using namespace sphlatch;

int main(int argc, char* argv[])
{
  const size_t noPartsX = 10;
  const size_t noPartsY = 10;
  const size_t noPartsZ = 100;

  const size_t noParts = noPartsX*noPartsY*noPartsZ;

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
  
  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.setNoParts(noParts);
  PartManager.resizeAll();

  std::cout << noParts << " particles \n";

  idvectRefType id( PartManager.id );
  
  matrixRefType pos( PartManager.pos );
  matrixRefType vel( PartManager.vel );
  
  valvectRefType m( PartManager.m );
  valvectRefType u( PartManager.u );
  valvectRefType h( PartManager.h );

  size_t curIndex = 0;

  for (size_t k = 0; k < noPartsZ; k++)
  {
    for (size_t j = 0; j < noPartsY; j++)
    {
      for (size_t i = 0; i < noPartsX; i++)
      {
        id(curIndex) = i + 10*j + 100*k;
        h(curIndex) = 1.;
        u(curIndex) = 1.;

        pos(curIndex, X) = static_cast<valueType>(i);
        pos(curIndex, Y) = static_cast<valueType>(j);
        pos(curIndex, Z) = static_cast<valueType>(k);

        vel(curIndex, X) = 0.;
        vel(curIndex, Y) = 0.;
        vel(curIndex, Z) = 0.;

        m(curIndex) = ( k > noPartsZ / 2 ) ? 1 : 2;
        
        curIndex++;
      }
    }
  }
  
  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &m, &u, &h;
  saveQuants.ints  += &id;

  PartManager.attributes["time"] = 0.0;
  PartManager.step = 0;
  
  //IOManager.setSinglePrecOut();
  IOManager.setDoublePrecOut();
  
  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cout << " -> " << outputFilename << "\n";
  IOManager.saveDump( outputFilename, saveQuants );

  return EXIT_SUCCESS;
}



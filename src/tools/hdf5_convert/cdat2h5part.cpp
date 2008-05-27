// some defs

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

#define OOSPH_STANDALONE
#define OOSPH_SINGLE_PRECISION

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

#include "simulation_trait.h"
typedef oosph::SimulationTrait<> SimTrait;
typedef SimTrait::matrix_reference oosph_matrix_ref_type;

#include "iomanager.h"
typedef oosph::IOManager<SimTrait> oosph_io_type;

#include "memorymanager.h"
typedef oosph::MemoryManager<SimTrait> oosph_mem_type;

#include "particle.h"

using namespace boost::assign;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i" , po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file") && !VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }
  
  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());
  oosph_io_type&  OOSPH_IOManager(oosph_io_type::Instance());
  oosph_mem_type& OOSPH_MemManager(oosph_mem_type::Instance());

  std::string outputFileName = VMap["output-file"].as<std::string>();
  std::string  inputFileName = VMap["input-file"].as<std::string>();

  OOSPH_IOManager.LoadCDAT(inputFileName);
  oosph_matrix_ref_type Data(OOSPH_MemManager.Data);

  const size_t noParts = Data.size1();

  PartManager.useGravity();
  PartManager.useBasicSPH();

  PartManager.setNoParts(noParts, 0);
  PartManager.resizeAll();

  sphlatch::rangeType partRange(0, noParts);

  sphlatch::matrixRefType pos( PartManager.pos );
  sphlatch::rangeType posRange(oosph::X, oosph::X+3);
  pos = sphlatch::matrixRangeType( Data, partRange, posRange );
  
  sphlatch::matrixRefType vel( PartManager.vel );
  sphlatch::rangeType velRange(oosph::VX, oosph::VX+3);
  vel = sphlatch::matrixRangeType( Data, partRange, posRange );

  sphlatch::valvectRefType h( PartManager.h );
  h   = sphlatch::quantColumnType( Data, oosph::H );
  
  sphlatch::valvectRefType m( PartManager.m );
  m   = sphlatch::quantColumnType( Data, oosph::M );
  
  sphlatch::valvectRefType eps( PartManager.eps );
  eps = sphlatch::quantColumnType( Data, oosph::GRAVEPS );

  sphlatch::idvectRefType id( PartManager.id );
  for (size_t i = 0; i < noParts; i++)
  {
    id(i) = lrint( Data(i, oosph::ID) );
  }

  using namespace boost::assign;
  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &h, &m, &eps;
  saveQuants.ints  += &id;

  PartManager.attributes["time"] = 0.0;
  PartManager.step = 0;
  
  //IOManager.setSinglePrecOut();
  IOManager.setDoublePrecOut();
  IOManager.saveDump( outputFileName, saveQuants );

  return EXIT_SUCCESS;
}



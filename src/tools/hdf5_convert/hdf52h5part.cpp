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

namespace po = boost::program_options;

#include "typedefs.h"
typedef sphlatch::valueType valueType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");
  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i" , po::value<std::string>(), "input  file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help") || not VMap.count("input-file") )
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  
  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  std::string  inputFilename = VMap["input-file"].as<std::string>();

  std::cout << "adding scalar position and velocity to " << inputFilename << "\n";

  IOManager.loadDump( inputFilename );

  //PartManager.useGravity();
  //PartManager.useBasicSPH();

  
  sphlatch::matrixRefType pos( PartManager.pos );
  sphlatch::matrixRefType vel( PartManager.vel );
  

  //IOManager.setSinglePrecOut();
  //IOManager.setDoublePrecOut();

  sphlatch::valvectType tmp;
  tmp.resize( pos.size1() );

  sphlatch::quantsType saveScalVects;
  saveScalVects.scalars += &tmp;

  PartManager.regQuantity( tmp, "pos_x" );
  tmp = sphlatch::quantColumnType( pos, 0 );
  IOManager.saveDump( inputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  PartManager.regQuantity( tmp, "pos_y" );
  tmp = sphlatch::quantColumnType( pos, 1 );
  IOManager.saveDump( inputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );

  PartManager.regQuantity( tmp, "pos_z" );
  tmp = sphlatch::quantColumnType( pos, 2 );
  IOManager.saveDump( inputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  PartManager.regQuantity( tmp, "vel_x" );
  tmp = sphlatch::quantColumnType( vel, 0 );
  IOManager.saveDump( inputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  PartManager.regQuantity( tmp, "vel_y" );
  tmp = sphlatch::quantColumnType( vel, 1 );
  IOManager.saveDump( inputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );

  PartManager.regQuantity( tmp, "vel_z" );
  tmp = sphlatch::quantColumnType( vel, 2 );
  IOManager.saveDump( inputFilename, saveScalVects );
  PartManager.unRegQuantity( tmp );
  
  return EXIT_SUCCESS;
}



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
typedef sphlatch::fType fType;
typedef sphlatch::stringListType stringListType;
typedef sphlatch::matrixRangeType matrixRangeType;
typedef sphlatch::rangeType rangeType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "simulation_trait.h"
typedef oosph::SimulationTrait<float> SimTrait;
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

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::string  inputFilename = VMap["input-file"].as<std::string>();

  stringListType loadVars = IOManager.discoverVars(inputFilename);

  stringListType::iterator varItr = loadVars.begin();
  stringListType::iterator varEnd = loadVars.end();

  using namespace oosph;
  size_t maxVar = Z;

  while( varItr != varEnd )
  {
    if ( *varItr == "rho")
    {
      PartManager.useBasicSPH();
      maxVar = RHO > maxVar ? RHO : maxVar;
    }
    
    if ( *varItr == "p")
    {
      PartManager.useBasicSPH();
      maxVar = P   > maxVar ? P   : maxVar;
    }
    
    if ( *varItr == "h")
    {
      PartManager.useBasicSPH();
      maxVar = H   > maxVar ? H   : maxVar;
    }
    
    if ( *varItr == "eps")
    {
      PartManager.useGravity();
      maxVar = GRAVEPS > maxVar ? GRAVEPS : maxVar;
    }
    
    if ( *varItr == "u")
    {
      PartManager.useEnergy();
      maxVar = E > maxVar ? E : maxVar;
    }
    
    if ( *varItr == "dudt")
    {
      PartManager.useTimedepEnergy();
      maxVar = POW > maxVar ? POW : maxVar;
    }
    
    if ( *varItr == "dhdt")
    {
      PartManager.useTimedepH();
      maxVar = DHDT > maxVar ? DHDT : maxVar;
    }
    
    if ( *varItr == "divv")
    {
      PartManager.useTimedepH();
      maxVar = DIVV > maxVar ? DIVV : maxVar;
    }
    
    if ( *varItr == "noneigh")
    {
      PartManager.useTimedepH();
      maxVar = NONEIGH > maxVar ? NONEIGH : maxVar;
    }
    
    if ( *varItr == "mumax")
    {
      PartManager.useAVMonaghan();
      maxVar = MUMAX > maxVar ? MUMAX : maxVar;
    }
    
    if ( *varItr == "alpha")
    {
      PartManager.useAVTimedepAlpha();
      maxVar = ALPHA > maxVar ? ALPHA : maxVar;
    }
    
    if ( *varItr == "dalphadt")
    {
      PartManager.useAVTimedepAlpha();
      maxVar = DALPHADT > maxVar ? DALPHADT : maxVar;
    }
    
    if ( *varItr == "rotv")
    {
      PartManager.useAVBalsara();
      maxVar = ROTZ_V > maxVar ? ROTZ_V : maxVar;
    }
    
    if ( *varItr == "q")
    {
      PartManager.useAVHernquist();
      maxVar = Q > maxVar ? Q : maxVar;
    }
    
    varItr++;
  }

  IOManager.loadDump(inputFilename);
  const size_t noParts = PartManager.getNoLocalParts();

  std::vector<int> saveVars;
  using namespace boost::assign;

  oosph_matrix_ref_type Data(OOSPH_MemManager.Data);
  sphlatch::rangeType partRange(0, noParts);

  std::cout << "DEBUG!!! " << PartManager.pos.size1() << "\n";

  varItr = loadVars.begin();
  while( varItr != varEnd )
  {
    if ( *varItr == "pos" )
    {
      std::cout << "pos ";
      rangeType posRange(X, X+3);
      matrixRangeType( Data, partRange, posRange ) = PartManager.pos;
      saveVars += X, Y, Z;
    }
    varItr++;
  }
  std::cout << "\n";
  
  OOSPH_IOManager.SaveCDAT(outputFilename, saveVars);

  /*sphlatch::matrixRefType vel( PartManager.vel );
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
  }*/

  /*using namespace boost::assign;
  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &h, &m, &eps;
  saveQuants.ints  += &id;

  PartManager.attributes["time"] = 0.0;
  PartManager.step = 0;
  
  //IOManager.setSinglePrecOut();
  IOManager.setDoublePrecOut();
  IOManager.saveDump( outputFileName, saveQuants );*/

  return EXIT_SUCCESS;
}



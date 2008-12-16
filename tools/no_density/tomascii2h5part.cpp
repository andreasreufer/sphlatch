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
typedef sphlatch::fType fType;

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

  if (!VMap.count("output-file") &&
      !VMap.count("input-file") )
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }
  
  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::string  inputFilename = VMap["input-file"].as<std::string>();

  std::cout << inputFilename << " -> " << outputFilename << "\n";

  std::fstream fin;
  fin.open(inputFilename.c_str(), std::ios::in);

  if (!fin)
  {
    std::cerr << "Error opening " << inputFilename << "\n";
  }

  size_t noParts = 0;
  fType gravConst = 0;
  std::string readLine;

  // ugly header parsing
  for (size_t i = 0; i < 8; i++)
  {
    fin >> readLine;
    std::cout << i << " " << readLine << "\n";
    switch (i)
    {
      case 1:
        noParts = boost::lexical_cast<size_t>(readLine);
        break;
    }
  }

  PartManager.useBasicSPH();
  PartManager.setNoParts(noParts, 0);
  PartManager.resizeAll();

  std::cout << noParts << " particles \n";

  sphlatch::idvectRefType id( PartManager.id );
  
  sphlatch::matrixRefType pos( PartManager.pos );
  
  sphlatch::valvectRefType h( PartManager.h );
  
  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos;
  saveQuants.scalars += &h;
  saveQuants.ints  += &id;

  PartManager.attributes["time"] = 0.0;
  PartManager.step = 0;
  
  for (size_t i = 0; i < noParts; i++)
  {
    for (size_t j = 0; j < 2; j++)
    {
      fin >> readLine;
      pos(i, j) = boost::lexical_cast<fType>(readLine);
    }

    pos(i, 2) = 0.;
   
    ///
    /// he's using search radius
    ///
    fin >> readLine;
    h(i) = 0.5*boost::lexical_cast<fType>(readLine);

    id(i) = i + 1;
  }

  fin.close();

  //IOManager.setSinglePrecOut();
  //IOManager.setDoublePrecOut();
  IOManager.saveDump( outputFilename, saveQuants );

  return EXIT_SUCCESS;
}



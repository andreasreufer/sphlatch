// some defs

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC
#define SPHLATCH_PARALLEL

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
typedef sphlatch::fType             fType;
typedef sphlatch::identType         identType;
typedef sphlatch::matrixType        matrixType;
typedef sphlatch::stringVectType    stringVectType;
typedef sphlatch::valvectType       valvectType;
typedef sphlatch::idvectType        idvectType;

#include "cdat_reader.h"
typedef sphlatch::CDATreader        cdatreader_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;


int main(int argc, char* argv[])
{
   MPI::Init(argc, argv);

   po::options_description Options("Global Options");
   Options.add_options()
   ("help,h", "Produces this Help")
   ("input-file,i", po::value<std::string>(), "input  file")
   ("output-file,o", po::value<std::string>(), "output file");

   po::variables_map VMap;
   po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
   po::notify(VMap);

   if (VMap.count("help"))
   {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
   }

   if (!VMap.count("output-file") && !VMap.count("input-file"))
   {
      std::cerr << Options << std::endl;
      return(EXIT_FAILURE);
   }

   cdatreader_type cdatReader;

   std::string outputFileName = VMap["output-file"].as<std::string>();
   std::string inputFilename  = VMap["input-file"].as<std::string>();

   cdatReader.open(inputFilename);

   const size_t noParts = cdatReader.getNoParts();
   const stringVectType vars = cdatReader.getVars();
   const size_t noVars = vars.size();

   std::fstream fout;
   fout.open(outputFileName.c_str(), std::ios::out);

   fout << "# " << noParts << " " << noVars << "\n";

   for (size_t j = 0; j < noVars; j++)
     fout << vars[j] << "\t";

   fout << "\n";

   for (size_t i = 0; i < noParts; i++)
   {
     for (size_t j = 0; j < noVars; j++)
     {
       fout << cdatReader.read(i,j) << "\t";
     }
     fout << "\n";
   }

   fout.close();

   cdatReader.close();

   MPI::Finalize();
   return(EXIT_SUCCESS);
}

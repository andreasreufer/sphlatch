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

using namespace boost::assign;
//using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file")
  ("header-file,f", po::value<std::string>(), "header file default: <none>");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help") ||
      not VMap.count("input-file") ||
      not VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }


  std::string inputFilename = VMap["input-file"].as<std::string>();
  std::string outputFilename = VMap["output-file"].as<std::string>();

  std::fstream fout;
  fout.open(outputFilename.c_str(), std::ios::out);

  if (!fout)
    {
      std::cerr << "Error Writing " << VMap[ "output-file" ].as<std::string>().c_str() << "\n";
    }

  if (VMap.count("header-file"))
    {
      std::string headerFilename = VMap["header-file"].as<std::string>();
      std::cout << " using header file " << headerFilename << "\n";

      std::string buff;
      std::fstream header;
      header.open(headerFilename.c_str(), std::ios::in);

      char c;
      while (header.get(c))
        {
          fout.put(c);
        }
      header.close();
    }
  else
    {
      std::cout << " using default header\n";
      fout << "#version 3.6;\n";
      fout << "global_settings {  assumed_gamma 1.0 }\n";
      fout << "#default{ finish{ ambient 0.15 diffuse 0.85}}\n";
      fout << "#include \"colors.inc\"\n";
      fout << "camera {location <3, 1, -5>\n";
      fout << "        look_at  < 0.5, 0 ,0>}\n";
      fout << "light_source{<15,0,-5> color White}\n";
    }

  std::fstream fin;
  fin.open(inputFilename.c_str(), std::ios::in);

  std::string curToken;

  const fType scale = 1.e-09;
  //const fType appRad = 0.002;  // value for 2Mparts
  const fType appRad = 0.010;    // value for 200kparts
  while (fin)
    {
      fout << "sphere{<";
      fin >> curToken;
      fout << scale*boost::lexical_cast<fType>(curToken);
      fout << ",";
      fin >> curToken;
      fout << scale*boost::lexical_cast<fType>(curToken);
      fout << ",";
      fin >> curToken;
      fout << scale*boost::lexical_cast<fType>(curToken);
      fout << ">, "<< appRad << "\n";
      fout << "   texture{\n";
      fout << "       pigment {color rgb<";
      fin >> curToken;
      fout << boost::lexical_cast<fType>(curToken);
      fout << ",";
      fin >> curToken;
      fout << boost::lexical_cast<fType>(curToken);
      fout << ",";
      fin >> curToken;
      fout << boost::lexical_cast<fType>(curToken);
      fout << ">} \n       } \n} \n";
    }

  fout.close();

  std::cout << inputFilename << " -> " << outputFilename << "\n";
  return EXIT_SUCCESS;
}




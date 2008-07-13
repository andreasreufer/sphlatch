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

  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file")
  ("header-file,f", po::value<std::string>(), "header file     default: <none>")
  ("scaling,s", po::value<valueType>(), "scaling         default: 1.00")
  ("radius,r", po::value<valueType>(), "sphere radius   default: 0.05");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help") || not VMap.count("input-file") || not VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  valueType k = 1.;
  if (VMap.count("scaling"))
    {
      k = VMap["scaling"].as<valueType>();
    }

  valueType rad = 0.05;
  if (VMap.count("radius"))
    {
      rad = VMap["radius"].as<valueType>();
    }
  rad *= k;

  std::cout << " scale particle positions by " << k
            << ", sphere radius is " << rad << "\n";

  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());

  std::string inputFilename = VMap["input-file"].as<std::string>();
  std::string outputFilename = VMap["output-file"].as<std::string>();

  PartManager.useBasicSPH();
  IOManager.loadDump(inputFilename);

  sphlatch::matrixRefType pos(PartManager.pos);
  sphlatch::valvectRefType h(PartManager.h);
  sphlatch::idvectRefType id(PartManager.id);

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
      fout << "#default{ finish{ ambient 0.1 diffuse 0.9}}\n";
      fout << "#include \"colors.inc\"\n";
      fout << "camera {location <0.7 , 0.7 ,-1.2>\n";
      fout << "        look_at  <0.0 , 0.0 , 0.0>}\n";
      fout << "light_source{<1500,2500,-2500> color White}\n";
    }

  const std::string redString = "rgb<1,0,0>";
  const std::string blueString = "rgb<0,0,1>";

  const size_t noParts = PartManager.getNoLocalParts();
  for (size_t i = 0; i < noParts; i++)
    {
      const valueType curX = k * pos(i, X);
      const valueType curY = k * pos(i, Y);
      const valueType curZ = k * pos(i, Z);

      const valueType sphereRad = 0.5 * k * h(i);

      if (curY < 0)
        {
          if (curY > -sphereRad)
            {
              fout << " difference {\n";
            }
          fout << "sphere{<"
               << curX << "," << curY << "," << curZ
               << ">, " << sphereRad << "\n";
          fout << "       pigment {color ";
          if (id(i) >= 1e6)
            fout << redString;
          else
            fout << blueString;
          fout << "}\n} \n";
          if (curY > -sphereRad)
            {
              const valueType boxSide =
                sqrt(sphereRad * sphereRad - curY * curY);
              fout << "box { <"
                   << curX - boxSide
                   << ",0,"
                   << curZ - boxSide
                   << "> , <"
                   << curX + boxSide
                   << ","
                   << boxSide
                   << ","
                   << curZ + boxSide
                   << "> pigment{ color ";
              if (id(i) >= 1e6)
                fout << redString;
              else
                fout << blueString;
              fout << " } }\n}\n";
            }
        }
    }

  fout.close();

  std::cout << inputFilename << " -> " << outputFilename << "\n";
  return EXIT_SUCCESS;
}




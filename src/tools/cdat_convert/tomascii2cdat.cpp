/**************************************************************************
 *   Copyright (C) 2006 by Andreas Reufer                                  *
 *   andreas.reufer@space.unibe.ch                                         *
 *                                                                         *
 ***************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef OOSPH_STANDALONE
#define OOSPH_STANDALONE
#endif

#include <iostream>
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

#include "simulation_trait.h"
#include "iomanager.h"
#include "particle.h"

namespace num = boost::numeric::ublas;
namespace po = boost::program_options;

typedef oosph::SimulationTrait<> SimTrait;
typedef oosph::IOManager<SimTrait> IOManagerType;
typedef oosph::MemoryManager<SimTrait> MemoryManagerType;
typedef oosph::ParticleVarMap<SimTrait> VarMapType;

typedef num::matrix<double>::size_type size_type;

typedef SimTrait::value_type value_type;
typedef SimTrait::vector_type vector_type;
typedef SimTrait::matrix_type matrix_type;
typedef SimTrait::matrix_reference matrix_reference;

int main(int argc, char* argv[])
{
  IOManagerType& IOManager(IOManagerType::Instance());
  MemoryManagerType& MemoryManager(MemoryManagerType::Instance());
  matrix_reference Data(MemoryManager.Data);

  po::options_description Options("Global Options");

  Options.add_options() ("help", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "Input   ASCII File")
  ("output-file,o", po::value<std::string>(), "Output  CDAT File")
  ("dataset-name,n", po::value<std::string>(), "Name of CDAT dataset         (default:   <none>)")
  ("grav-epsilon,e", po::value<value_type>(), "gravitational smoothing      (default:       0 )");

  po::options_description Config("Config File Options");

  po::positional_options_description p;
  p.add("input-file", 1);
  p.add("output-file", 2);
  p.add("dataset-name", 3);
  p.add("grav-epsilon", 4);

  po::variables_map VMap;

  po::store(po::command_line_parser(argc, argv).options(Options).positional(p).run(), VMap);

  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cout << Options << std::endl;
      exit(-1);
    }

  if (!VMap.count("input-file"))
    {
      std::cout << std::endl << "Input File not Specified" << std::endl << std::endl;
      std::cout << Options << std::endl;
      exit(-1);
    }
  std::string inputFilename = VMap[ "input-file" ].as<std::string>();

  if (!VMap.count("output-file"))
    {
      std::cout << std::endl << "Output File has to be specified" << std::endl << std::endl;
      std::cout << Options << std::endl;
      exit(-1);
    }
  std::string outputFilename = VMap[ "output-file" ].as<std::string>();

  if (!VMap.count("dataset-name"))
    {
      MemoryManager.Name = "<none>";
    }
  else
    {
      MemoryManager.Name = VMap[ "dataset-name" ].as<std::string>();
    }

  value_type gravEps = 0.;
  if (VMap.count("grav-epsilon"))
    {
      gravEps = VMap[ "grav-epsilon" ].as<value_type>();
    }

  // **
  // ** Start Reading in Data
  // **

  std::cout << inputFilename << " -> " << outputFilename << "\n";

  std::fstream fin;
  fin.open(inputFilename.c_str(), std::ios::in);

  if (!fin)
    {
      std::cerr << "Error opening " << inputFilename << "\n";
    }

  using boost::lexical_cast;
  using namespace boost::assign;
  using namespace oosph;

  std::string readLine;

  size_t noParts;
  value_type gravConst;

  // ugly header parsing
  for (size_t i = 0; i < 25; i++)
    {
      fin >> readLine;
      switch (i)
        {
        case 1:
          noParts = lexical_cast<size_t>(readLine);
          break;

        case 5:
          gravConst = lexical_cast<value_type>(readLine);
          break;
        }
    }

  std::cout << noParts << " particles \n";
  std::cout << "using G = " << gravConst
            << " and gravitational smoothing of "
            << gravEps << "\n";

  Data.resize(noParts, GRAVEPS + 1);

  const size_t noVars = 7;    // m, x, y, z, vx, vy, vz
  std::vector<int> outputAttrSet;
  outputAttrSet += M, X, Y, Z, VX, VY, VZ;

  for (size_t i = 0; i < noParts; i++)
    {
      for (size_t j = 0; j < noVars; j++)
        {
          fin >> readLine;
          Data(i, outputAttrSet[j]) = lexical_cast<value_type>(readLine);
        }
      Data(i, ID) = lexical_cast<value_type>(i + 1);
      Data(i, GRAVEPS) = gravEps;
    }

  MemoryManager.SaveParameter("TIME", 0, true);
  MemoryManager.SaveParameter("GRAVCONST", gravConst, true);

  fin.close();

  outputAttrSet += ID, GRAVEPS;
  IOManager.SaveCDAT(outputFilename, outputAttrSet);

  std::cout << "saved to " <<  outputFilename << "!\n";

  exit(EXIT_SUCCESS);
}


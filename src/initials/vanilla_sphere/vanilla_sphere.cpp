// some defs

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
typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "lattice_hcp.h"
typedef sphlatch::LatticeHCP body_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options()
  ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), "output file")
  ("no-parts,n",    po::value<size_t>(),      "number of particles");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file") || !VMap.count("no-parts"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type&        IOManager(io_type::instance());
  part_type&      PartManager(part_type::instance());


  matrixRefType pos(PartManager.pos);
  valvectRefType  m(PartManager.m);
  valvectRefType  h(PartManager.h);
  idvectRefType  id(PartManager.id);

  const size_t desNoParts = VMap["no-parts"].as<size_t>();

  ///
  /// 1.03 is a fudge factor :-)
  /// 
  const valueType rMax = 1.;
  const valueType fillingFactor = 0.740;
  const valueType latticeLength = 2.*rMax/pow( fillingFactor*desNoParts
                                               / 1.03, 1./3.);

  ///
  /// now place the SPH particles on a lattice
  ///
  sphlatch::LatticeHCP Lattice(latticeLength, 1.1*rMax, 1.1*rMax, 1.1*rMax);

  size_t partsCount = 0;
  Lattice.first();
  while (!Lattice.isLast)
    {
      if (Lattice.rCur < rMax)
        partsCount++;
      Lattice.next();
    }
  std::cerr << "you asked for " << desNoParts
            << " and you get " << partsCount << " particles\n";

  PartManager.useBasicSPH();
  PartManager.setNoParts(partsCount);
  PartManager.resizeAll();
  
  const valueType unitVolume = pow(latticeLength,3.)*fillingFactor;
  const valueType smoLength  = 1.*latticeLength;

  std::cout << "tot volume " << partsCount*unitVolume << "\n";

  partsCount = 0;
  Lattice.first();
  while (!Lattice.isLast)
    {
      if (Lattice.rCur < rMax)
        {
          id(partsCount) = partsCount;

          pos(partsCount, X) = Lattice.xCur;
          pos(partsCount, Y) = Lattice.yCur;
          pos(partsCount, Z) = Lattice.zCur;

          ///
          /// the unit volume is stored in the mass variable
          ///
          /// so, multiplied with the density, one gets the
          /// particle mass
          ///
          m(partsCount) = unitVolume;
          h(partsCount) = smoLength;

          partsCount++;
        }
      Lattice.next();
    }

  sphlatch::quantsType saveQuants;
  saveQuants.vects += &pos;
  saveQuants.scalars += &m, &h;
  saveQuants.ints += &id;
  PartManager.step = 0;

  std::string outputFilename = VMap["output-file"].as<std::string>();
  std::cerr << " -> " << outputFilename << "\n";
  IOManager.saveDump(outputFilename, saveQuants);

  return EXIT_SUCCESS;
}



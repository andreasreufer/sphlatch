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
typedef sphlatch::fType fType;
typedef sphlatch::valvectType valvectType;

typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::quantsType quantsType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

extern "C" void writecracks_(float *epsmin, float *xm,
                             int *nflaws, float *acoef,
                             float *young, float *grav);

#include "cdat_writer.h"
typedef sphlatch::CDATwriter cdatwriter_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");

  Options.add_options()
  ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help") ||
      not VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }


  part_type& PartManager(part_type::instance());
  io_type&        IOManager(io_type::instance());
  cdatwriter_type CDATWriter;

  std::string inputFilename = VMap["input-file"].as<std::string>();

  std::cerr << inputFilename << "\n"
            << " -> initial.xdr\n";

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useMaterials();
  PartManager.useDamage();
  PartManager.useStress();

  IOManager.loadDump(inputFilename);

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  matrixRefType S(PartManager.S);

  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType rho(PartManager.rho);
  valvectRefType u(PartManager.u);
  valvectRefType p(PartManager.p);
  valvectRefType dam(PartManager.dam);
  valvectRefType epsmin(PartManager.epsmin);
  valvectRefType acoef(PartManager.acoef);
  valvectRefType mweib(PartManager.mweib);
  valvectRefType young(PartManager.young);

  idvectRefType id(PartManager.id);
  idvectRefType mat(PartManager.mat);
  idvectRefType noflaws(PartManager.noflaws);

  const size_t noParts = PartManager.getNoLocalParts();

  std::vector<std::string> varNames;

  varNames += "x", "y", "z", "vx", "vy", "vz", "h", "id", "mass",
  "mat", "p", "T", "rho", "u", "dm", "rft", "sxx", "sxy",
  "sxz", "syy", "syz";

  CDATWriter.open("initial.xdr", "none", PartManager.attributes["time"],
                  noParts, varNames, false);

  for (size_t i = 0; i < noParts; i++)
    {
      CDATWriter.write(i, 0, pos(i, X));
      CDATWriter.write(i, 1, pos(i, Y));
      CDATWriter.write(i, 2, pos(i, Z));

      CDATWriter.write(i, 3, vel(i, X));
      CDATWriter.write(i, 4, vel(i, Y));
      CDATWriter.write(i, 5, vel(i, Z));

      CDATWriter.write(i, 6, h(i));
      CDATWriter.write(i, 7, static_cast<fType>(id(i)));
      CDATWriter.write(i, 8, m(i));
      CDATWriter.write(i, 9, static_cast<fType>(mat(i)));
      CDATWriter.write(i, 10, p(i));
      CDATWriter.write(i, 11, 0.);
      CDATWriter.write(i, 12, rho(i));
      CDATWriter.write(i, 13, u(i));
      CDATWriter.write(i, 14, dam(i));
      CDATWriter.write(i, 15, 0.);

      CDATWriter.write(i, 16, S(i, XX));
      CDATWriter.write(i, 17, S(i, XY));
      CDATWriter.write(i, 18, S(i, XZ));
      CDATWriter.write(i, 19, S(i, YY));
      CDATWriter.write(i, 20, S(i, YZ));
    }

  CDATWriter.close();

  std::cerr << " -> in.cracks\n";
  for (size_t i = 0; i < noParts; i++)
    {
      float curEpsmin = epsmin(i);
      float curMweib = mweib(i);
      int curNoflaws = noflaws(i);
      float curAcoef = acoef(i);
      float curYoung = young(i);
      float curGravsh = 0.;

      writecracks_(&curEpsmin, &curMweib, &curNoflaws,
                   &curAcoef, &curYoung, &curGravsh);
    }

  return EXIT_SUCCESS;
}


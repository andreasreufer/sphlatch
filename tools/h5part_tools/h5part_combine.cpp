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
typedef sphlatch::fType fType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);
  po::options_description Options("Global Options");

  Options.add_options() ("help,h", "Produces this Help")
  ("inputa-file,a", po::value<std::string>(), "input file A")
  ("inputb-file,b", po::value<std::string>(), "inout file B")
  ("output-file,o", po::value<std::string>(), "output file");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (!VMap.count("output-file") || VMap.count("help") ||
      !VMap.count("inputa-file") || !VMap.count("inputb-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  std::string inputAFileName = VMap["inputa-file"].as<std::string>();
  std::string inputBFileName = VMap["inputb-file"].as<std::string>();
  std::string outputFileName = VMap["output-file"].as<std::string>();

  PartManager.useBasicSPH();
  PartManager.useMaterials();
  PartManager.useEnergy();
  PartManager.useGravity();
#ifdef SPHLATCH_SOLID
  PartManager.useDamage();
  PartManager.useStress();
#endif

  using namespace sphlatch;
  using namespace boost::assign;

  idvectRefType id(PartManager.id);
  idvectRefType mat(PartManager.mat);
#ifdef SPHLATCH_SOLID
  idvectRefType noflaws(PartManager.noflaws);
#endif

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
#ifdef SPHLATCH_SOLID
  matrixRefType S(PartManager.S);
#endif

  valvectRefType m(PartManager.m);
  valvectRefType u(PartManager.u);
  valvectRefType h(PartManager.h);
  valvectRefType eps(PartManager.eps);
#ifdef SPHLATCH_SOLID
  valvectRefType rho(PartManager.rho);
  valvectRefType dam(PartManager.dam);
  valvectRefType epsmin(PartManager.epsmin);
  valvectRefType acoef(PartManager.acoef);
  valvectRefType mweib(PartManager.mweib);
  valvectRefType young(PartManager.young);
#endif

  quantsType saveQuants;
  saveQuants.ints += &id, &mat;
  saveQuants.vects += &pos, &vel;
  saveQuants.scalars += &m, &u, &h, &eps;
#ifdef SPHLATCH_SOLID
  saveQuants.ints += &noflaws;
  saveQuants.vects += &S;
  saveQuants.scalars += &dam, &epsmin, &acoef, &mweib, &young, &rho;
#endif

  IOManager.loadDump(inputAFileName);
  IOManager.loadDump(inputBFileName);
  IOManager.saveDump(outputFileName, saveQuants);

  MPI::Finalize();
  return EXIT_SUCCESS;
}



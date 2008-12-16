#define SPHLATCH_PARALLEL

#include <iostream>

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
typedef sphlatch::identType identType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::valvectType valvectType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::quantsType quantsType;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "parse_po_vectors.h"

using namespace boost::assign;
using namespace sphlatch::vectindices;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  po::options_description Options("Global Options");

  Options.add_options() ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input file")
  ("id,d", po::value<identType>(), "displace IDs")
  ("pos", po::value<std::string>(), "displace position ( += [x,y,z])")
  ("vel", po::value<std::string>(), "displace velocity ( += [vx,vy,vz])");

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), VMap);
  po::notify(VMap);

  if (!VMap.count("input-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());

  std::string inputFileName = VMap["input-file"].as<std::string>();
  quantsType saveQuants;

  PartManager.useBasicSPH();
  PartManager.useEnergy();
  PartManager.useGravity();
#ifdef SPHLATCH_SOLID
  PartManager.useMaterials();
  PartManager.useDamage();
  PartManager.useStress();
#endif

  IOManager.loadDump(inputFileName);
  const size_t noParts = PartManager.getNoLocalParts();

  bool changed = false;

  if (VMap.count("pos"))
    {
      matrixRefType pos(PartManager.pos);
      valvectType dpos = vectOptParse(VMap[ "pos" ].as<std::string>());

      std::cerr << "displacing position by ["
                << dpos(X) << ","
                << dpos(Y) << ","
                << dpos(Z) << "]\n";

      for (size_t i = 0; i < noParts; i++)
        {
          pos(i, X) += dpos(X);
          pos(i, Y) += dpos(Y);
          pos(i, Z) += dpos(Z);
        }

      saveQuants.vects += &pos;
      changed = true;
    }

  if (VMap.count("vel"))
    {
      matrixRefType vel(PartManager.vel);
      valvectType dvel = vectOptParse(VMap[ "vel" ].as<std::string>());

      std::cerr << "displacing velocity by ["
                << dvel(X) << ","
                << dvel(Y) << ","
                << dvel(Z) << "]\n";

      for (size_t i = 0; i < noParts; i++)
        {
          vel(i, X) += dvel(X);
          vel(i, Y) += dvel(Y);
          vel(i, Z) += dvel(Z);
        }

      saveQuants.vects += &vel;
      changed = true;
    }

  if (VMap.count("id"))
    {
      idvectRefType id(PartManager.id);
      identType did = VMap[ "id" ].as<identType>();

      std::cerr << "displacing ID by       " << did << "\n";

      for (size_t i = 0; i < noParts; i++)
        {
          id(i) += did;
        }

      saveQuants.ints += &id;
      changed = true;
    }

  if (changed)
    IOManager.saveDump(inputFileName, saveQuants);
  
  MPI::Finalize();
  return EXIT_SUCCESS;
}



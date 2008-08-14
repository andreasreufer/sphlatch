// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
namespace po = boost::program_options;

#include <boost/mpl/vector_c.hpp>

#include <boost/assign/std/vector.hpp>

#include "typedefs.h"
typedef sphlatch::valueType valueType;
typedef sphlatch::zerovalvectType zerovalvectType;
typedef sphlatch::valvectRefType valvectRefType;
typedef sphlatch::idvectRefType idvectRefType;
typedef sphlatch::matrixRefType matrixRefType;
typedef sphlatch::partsIndexVectType partsIndexVectType;
typedef sphlatch::quantsType quantsType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager comm_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include <boost/progress.hpp>
#include <vector>

// tree stuff
#include "bhtree.h"

using namespace sphlatch::vectindices;
using namespace boost::assign;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  po::options_description Options("Global Options");
  Options.add_options() ("help", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "hdf5 file")
  ("no-neighs,n", po::value<size_t>(), "no. of neighbours  (default 50)")
  ("cs-length,c", po::value<valueType>(), "comp. supp. length (default 2.)");

  po::variables_map PoMap;
  po::store(po::command_line_parser(argc, argv).options(Options).run(), PoMap);
  po::notify(PoMap);

  if (not PoMap.count("input-file"))
    {
      std::cerr << "input-file not specified!\n";
      std::cerr << Options << "\n";
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  std::string filename = PoMap[ "input-file" ].as<std::string>();

  size_t noDesNeighs = 50;
  if (PoMap.count("no-neighs"))
    {
      noDesNeighs = PoMap[ "no-neighs" ].as<size_t>();
    }

  valueType csLength = 2.;
  if (PoMap.count("cs-length"))
    {
      csLength = PoMap[ "cs-length" ].as<valueType>();
    }

  std::cout << " trying to get " << noDesNeighs << " neighbours inside "
            << csLength << "h\n";

  ///
  /// instantate managers
  ///
  io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  comm_type& CommManager(comm_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  ///
  /// define what we're doing
  ///
  PartManager.useBasicSPH();

  ///
  /// some useful references
  ///
  matrixRefType pos(PartManager.pos);
  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  ///
  /// register the quantites to be exchanged
  ///
  CommManager.exchangeQuants.vects += &pos;
  CommManager.exchangeQuants.scalars += &m;
  CommManager.exchangeQuants.ints += &id;

  IOManager.loadDump(filename);

  ///
  /// exchange particles
  ///
  CostZone.createDomainPartsIndex();
  CommManager.exchange(CostZone.domainPartsIndex,
                       CostZone.getNoGhosts());

  ///
  /// prepare ghost sends and send ghosts
  ///
  CommManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  CommManager.sendGhosts(pos);
  CommManager.sendGhosts(id);

  ///
  /// build empty tree
  ///
  sphlatch::BHtree<sphlatch::Monopoles>   Tree(1.0,
                                               0.0,
                                               CostZone.getDepth(),
                                               CostZone.getCenter(),
                                               CostZone.getSidelength()
                                               );


  ///
  /// fill up tree, determine ordering, calculate multipoles
  ///
  const size_t noParts = PartManager.getNoLocalParts();
  const size_t noTotParts = PartManager.getNoTotalParts();
  for (size_t i = 0; i < noTotParts; i++)
    {
      m(i) = 1.;
      Tree.insertParticle(i);
    }
  Tree.detParticleOrder();
  Tree.calcMultipoles();

  ///
  /// this is the maximal number of neighbours to be allowed
  /// in the maxmimal mass enclosing sphere. if the number's
  /// bigger, the neighbour search algorithm will unspecta-
  /// cularly crash
  ///
  const size_t noMaxNeighs = 16834;

  Tree.neighbourList.resize(noMaxNeighs);
  Tree.neighDistList.resize(noMaxNeighs);

  zerovalvectType zeroVect(noMaxNeighs);

  const valueType minMass = static_cast<valueType>(noDesNeighs);

  for (size_t i = 0; i < noParts; i++)
    {
      const size_t curIndex = Tree.particleOrder[i];

      const valueType maxRad = Tree.maxMassEncloseRad(curIndex, minMass);

      ///
      /// the tree does not guarantee valid entries in the neighbour list
      /// and the neighbour dist list above the number of actual neighbours,
      /// so zero the dist list before searching the neighbours
      ///
      Tree.neighDistList = zeroVect;
      Tree.findNeighbours(curIndex, maxRad);

      const size_t noGuessNeighs = Tree.neighbourList[0];

      ///
      /// sort neighbour list in ascending order
      ///
      std::sort(Tree.neighDistList.begin(), Tree.neighDistList.end());

      assert(noMaxNeighs - noGuessNeighs + noDesNeighs < noMaxNeighs);
      h(curIndex) = 1.01*Tree.neighDistList(noMaxNeighs - 1 - noGuessNeighs
                                       + noDesNeighs) / csLength;
      
      Tree.findNeighbours(curIndex, csLength*h(curIndex));
      noneigh(curIndex) = Tree.neighbourList[0];
    }

  ///
  /// overwrite smoothing length in file
  ///
  quantsType saveQuants;
  saveQuants.scalars += &h;
  saveQuants.ints += &noneigh;
  IOManager.saveDump(filename, saveQuants);

  MPI::Finalize();
  return EXIT_SUCCESS;
}



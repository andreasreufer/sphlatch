// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
//#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define GRAVITY
#define TREEORDER
#define RSORDER

//#define SPHLATCH_RANKSPACESERIALIZE

//#define BFCHECK
//#define CHECK_TREE
//#define CHECK_RANKSPACE

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

#include <boost/assign/std/vector.hpp>

#include <boost/mpl/vector_c.hpp>

namespace po = boost::program_options;

#include "typedefs.h"
//typedef sphlatch::fType fType;
//typedef sphlatch::partsIndexVectType partsIndexVectType;

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager com_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

#include <boost/progress.hpp>
#include <vector>

// tree stuff
#include "bhtree.h"

//#include "ranklist.h"
#include "rankspace.h"

#include "timer.h"

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);
  double stepStartTime, logStartTime = MPI_Wtime();

  po::options_description Options("Global Options");
  Options.add_options() ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input  file")
  ("output-file,o", po::value<std::string>(), "output file");

  po::positional_options_description POD;
  POD.add("input-file", 1);

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(POD).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  if (!VMap.count("input-file") or
      !VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  using namespace boost::assign;
  using namespace sphlatch;

  io_type& IOManager(io_type::instance());
  part_type& PartManager(part_type::instance());
  com_type& ComManager(com_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  const size_t myDomain = ComManager.getMyDomain();

  size_t noParts, noGhosts, noTotParts;

  std::string inputFileName = VMap[ "input-file"].as<std::string>();
  std::string outputFileName = VMap["output-file"].as<std::string>();

  PartManager.useBasicSPH();
  PartManager.useTimedepH();
  IOManager.loadDump(inputFileName);

  // set up logging stuff
  std::string logFilename = outputFileName + "_rank000";
  std::string rankString = boost::lexical_cast<std::string>(myDomain);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "    start log\n";

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);

  valvectRefType h(PartManager.h);
  valvectRefType m(PartManager.m);

  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  quantsType exchQuants;
  exchQuants.vects += &pos, &vel;
  exchQuants.ints += &id, &noneigh;
  exchQuants.scalars += &m, &h;

  stepStartTime = MPI_Wtime();
  /// domain decomposition and exchange of particles
  CostZone.createDomainPartsIndex();
  ComManager.exchange(CostZone.domainPartsIndex,
                      CostZone.getNoGhosts(),
                      exchQuants);

  /// prepare ghost stuff
  ComManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());
  ComManager.sendGhosts(exchQuants);

  noParts = PartManager.getNoLocalParts();
  noTotParts = PartManager.getNoTotalParts();

  sleep(5);

  size_t curIndex = 0;
  double treesearchTime, rssearchTime;

#ifdef RSORDER
  partsIndexVectType orderFromTree;
  orderFromTree.resize(noParts);
#endif

  // tree context
  {
    const fType gravTheta = 0.6;
    const fType gravConst = 1.0;

    stepStartTime = MPI_Wtime();
    
    //sphlatch::BHtree<sphlatch::NoMultipoles> BarnesHutTree(gravTheta,
    //sphlatch::BHtree<sphlatch::Monopoles> BarnesHutTree(gravTheta,
    sphlatch::BHtree<sphlatch::Quadrupoles> BarnesHutTree(gravTheta,
    //sphlatch::BHtree<sphlatch::Octupoles> BarnesHutTree(gravTheta,
                                                          gravConst,
                                                          //CostZone.getDepth(),
                                                          0,
                                                          CostZone.getCenter(),
                                                          CostZone.getSidelength());
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - stepStartTime << "lps constructed tree\n" << std::flush;

    for (size_t i = 0; i < noTotParts; i++)
    {
      BarnesHutTree.insertParticle(i);
    }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - stepStartTime << "lps inserted parts\n" << std::flush;


#ifdef TREEORDER
    BarnesHutTree.detParticleOrder();
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - stepStartTime << "lps det. part. order\n" << std::flush;
#endif

    for (size_t i = 0; i < noParts; i++)
      {
#ifdef TREEORDER
        curIndex = BarnesHutTree.particleOrder[i];
#else
        curIndex = i;
#endif
        // allow a little tolerance, as the furthest is exactely at 2h
        const fType searchRadius = 2.0000001 * h(curIndex);
        BarnesHutTree.findNeighbours(curIndex, searchRadius);
        std::cout << curIndex << ": " << BarnesHutTree.neighbourList[0] << "\n";
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - stepStartTime << "    BHtree search       \n" << std::flush;
    treesearchTime = MPI_Wtime() - stepStartTime;
    logFile << "\n";

#ifdef RSORDER
    orderFromTree = BarnesHutTree.particleOrder;
#endif

  }
  // end of tree context

  sleep(5);

  std::cout << "\n\n";

  // rankspace context
  {
    stepStartTime = MPI_Wtime();
    sphlatch::Rankspace RSSearch;
    RSSearch.prepare();
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - stepStartTime << "lps prepare RSSearch    \n" << std::flush;
    for (size_t i = 0; i < noParts; i++)
      {
#ifdef RSORDER
        curIndex = orderFromTree[i];
#else
        curIndex = i;
#endif
        // allow a little tolerance, as the furthest is exactely at 2h
        const fType searchRadius = 2.0000001 * h(curIndex);
        RSSearch.findNeighbours(curIndex, searchRadius);
        std::cout << curIndex << ": " << RSSearch.neighbourList[0] << "\n";
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
            << MPI_Wtime() - stepStartTime << "    Ranksspace search   \n" << std::flush;
    rssearchTime = MPI_Wtime() - stepStartTime;
  }
  // end of rankspace context

  logFile << "\n" << "  RS " << treesearchTime / rssearchTime << "x faster than tree \n";

#ifdef GRAVITY
  stepStartTime = MPI_Wtime();
  BarnesHutTree.calcMultipoles();
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime << "    calculated multipoles\n" << std::flush;

  for (size_t i = 0; i < noParts; i++)
    {
      curIndex = BarnesHutTree.particleOrder[i];
      //curIndex = i;
      //BarnesHutTree.calcGravity(*(partProxies[curIndex]));
      BarnesHutTree.calcGravity(i);
    }
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime << "    calculated gravity\n" << std::flush;
#endif

#ifdef BFCHECK
  partsIndexVectType neighListBF;
  neighListBF.resize(100);

  stepStartTime = MPI_Wtime();
  for (size_t i = 0; i < noParts; i++)
  //for (size_t i = 0; i < 100; i++)
  //size_t i = 88358;
    {
      curIndex = i;

      // allow a little tolerance, as the furthest is exactely at 2h
      const fType searchRadius = 2.0000001 * h(curIndex);
      const fType searchRadPow2 = searchRadius * searchRadius;

      size_t noBFneighbours = 0;
      for (size_t j = 0; j < noParts; j++)
        {
          static fType dist;

          //using namespace sphlatch;
          /*dist = sqrt((Data(i, X) - Data(j, X)) * (Data(i, X) - Data(j, X))
                      + (Data(i, Y) - Data(j, Y)) * (Data(i, Y) - Data(j, Y))
                      + (Data(i, Z) - Data(j, Z)) * (Data(i, Z) - Data(j, Z)));*/
          dist = sqrt((pow(i, X) - pow(j, X)) * (pow(i, X) - pow(j, X))
                      + (pow(i, Y) - pow(j, Y)) * (pow(i, Y) - pow(j, Y))
                      + (pow(i, Z) - pow(j, Z)) * (pow(i, Z) - pow(j, Z)));
          if (dist < searchRadius)
            {
              noBFneighbours++;
              neighListBF[noBFneighbours] = j;
            }
        }
      neighListBF[0] = noBFneighbours;

#ifdef CHECK_TREE
      //BarnesHutTree.findNeighbours(*(partProxies[curIndex]), searchRadius);
      BarnesHutTree.findNeighbours(curIndex, searchRadius);
      sphlatch::partsIndexVectRefType neighList(BarnesHutTree.neighbourList);
#endif
#ifdef CHECK_RANKSPACE
      RSSearch.findNeighbours(curIndex, searchRadius);
      sphlatch::partsIndexVectRefType neighList(RSSearch.neighbourList);
#endif
      size_t noBFneighs = neighListBF[0];
      neighListBF[0] = 999999;
      std::sort(neighListBF.begin(), neighListBF.end(), std::greater<size_t>());
      size_t noNeighbours = neighList[0];
      neighList[0] = 999999;
      std::sort(neighList.begin(), neighList.end(), std::greater<size_t>());

      bool neighboursEqual = true;
      for (size_t j = 0; j < 100; j++)
        {
          if (neighListBF[j] != neighList[j])
            {
              neighboursEqual = false;
              break;
            }
        }

      if (!neighboursEqual)
        {
          //std::cout << i << ":  " << noBFneighs << " (BF)  vs.  " << noNeighbours << " (tree)\n";
          std::cout << " tree search and BF search MISMATCH!!! (" << i << ")\n";
          std::cout << " BF    vs.     tree \n";
          for (size_t j = 0; j < 100; j++)
            {
              if (neighListBF[j] != neighList[j])
                {
                  std::cout << j << ":  " << neighListBF[j] << "      " << neighList[j] << "\n";
                }
            }
          std::cout << "\n";
        }
      else
        {
          std::cout << i << "\n";
        }
    }

  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime << " BF neighbour search\n" << std::flush;
#endif
  
  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}



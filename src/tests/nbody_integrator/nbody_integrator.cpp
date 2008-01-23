// some defs

#define OOSPH_SINGLE_PRECISION
#define SPHLATCH_SINGLEPREC

//#define OOSPH_LOADBALANCE
#define SPHLATCH_MPI

#ifdef SPHLATCH_MPI
#include <mpi.h>
#endif

#include <cstdlib>
#include <iostream>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include <boost/assign/std/vector.hpp>

#include <boost/mpl/vector_c.hpp>

#include "particle.h"

namespace po = boost::program_options;
namespace mpl = boost::mpl;

#include "simulation_trait.h"
typedef oosph::SimulationTrait<> SimTrait;
typedef SimTrait::value_type value_type;

#include "iomanager.h"
typedef oosph::IOManager<SimTrait> io_type;

#include "memorymanager.h"
typedef oosph::MemoryManager<SimTrait> mem_type;

#include "communicationmanager.h"
typedef oosph::CommunicationManager<SimTrait> com_type;

#include "costzone.h"
typedef mpl::vector_c<size_t, oosph::X> CostZoneIndex;
typedef oosph::CostZone<CostZoneIndex, SimTrait> CostZoneType;

#include "predictorcorrector.h"
typedef mpl::vector_c<size_t, oosph::X, oosph::VX, oosph::AX, oosph::OX, oosph::OVX, oosph::OAX, oosph::PX, oosph::PVX, oosph::PAX> AccIntIndices;
typedef oosph::SOVecPredictorCorrector<AccIntIndices, SimTrait> AccIntType;

#include "erazer.h"
typedef mpl::vector_c<size_t, oosph::AX, oosph::AY, oosph::AZ> ErazerIndices;
typedef oosph::Erazer<ErazerIndices, SimTrait> ErazerType;

using namespace oosph;
using namespace boost::assign;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <vector>

// tree stuff
#include "bhtree.h"

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  //const size_t SIZE = MPI::COMM_WORLD.Get_size();

  po::options_description Options("Global Options");
  Options.add_options() ("help,h", "Produces this Help blabla...")
  ("input-file,i", po::value<std::string>(), "InputFile");

  po::positional_options_description POD;
  POD.add("input-file", 1);

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(POD).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cout << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  if (!VMap.count("input-file"))
    {
      std::cout << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  io_type& IOManager(io_type::Instance());
  mem_type& MemManager(mem_type::Instance());
  com_type& ComManager(com_type::Instance());
  CostZoneType& CostZone(CostZoneType::Instance());

  ErazerType& Erazer(ErazerType::Instance());
  AccIntType AccInt;

  using namespace boost::posix_time;
  ptime TimeStart, TimeStop;

  std::vector<sphlatch::NodeProxy> partProxies, ghostProxies;

  SimTrait::matrix_reference Data(MemManager.Data);
  SimTrait::matrix_reference GData(MemManager.GData);

  Data.resize(Data.size1(), oosph::SIZE);
  GData.resize(GData.size1(), oosph::OX);

  std::string InputFileName = VMap["input-file"].as<std::string>();
  IOManager.LoadCDAT(InputFileName);
  value_type dt, absTime = 0.; // load from file
  value_type theta = 0.7;
  size_t noParts, noGhosts, step = 0;

  std::string logFilename = "logRank000";
  std::string rankString = boost::lexical_cast<std::string>(RANK);
  logFilename.replace(logFilename.size() - 3 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << MPI_Wtime() << "    start log\n";
  
  // bootstrapping context
  if (RANK == 0)
    {
      std::cerr << "-------- bootstrapping step -----\n";
    }
  {
    logFile << MPI_Wtime() << "    bootstrap\n";
    ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
    logFile << MPI_Wtime() << "    exchange parts\n";
    ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
    logFile << MPI_Wtime() << "    exchange ghosts\n";
    noParts = Data.size1();
    noGhosts = GData.size1();
    std::cerr << "Rank " << RANK << ": "
              << noParts << " particles and " << noGhosts << " ghosts\n";

    partProxies.resize(noParts);
    for (size_t i = 0; i < noParts; i++)
      {
        Data(i, M) = Data(i, M) / 100.; // adapt mass
        (partProxies[i]).setup(&Data, i);
      }
    ghostProxies.resize(noParts);
    for (size_t i = 0; i < noGhosts; i++)
      {
        GData(i, M) = GData(i, M) / 100.; // adapt mass
        (ghostProxies[i]).setup(&GData, i);
      }

    TimeStart = microsec_clock::local_time();
    sphlatch::OctTree BarnesHutTree(theta,
                                    CostZone.getDepth(),
                                    CostZone.getCenter(),
                                    CostZone.getSidelength());
    for (size_t i = 0; i < noParts; i++)
      {
        BarnesHutTree.insertParticle(*(partProxies[i]), true);
      }
    for (size_t i = 0; i < noGhosts; i++)
      {
        BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
      }
    BarnesHutTree.calcMultipoles();
    TimeStop = microsec_clock::local_time();
    std::cerr << "Rank " << RANK << ": "
              << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

    TimeStart = microsec_clock::local_time();
    for (size_t i = 0; i < noParts; i++)
      {
        BarnesHutTree.calcGravity(*(partProxies[i]));
      }
    TimeStop = microsec_clock::local_time();
    std::cerr << "Rank " << RANK << ": "
              << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
  }
  AccInt.BootStrap();
  Erazer();

  while (step < 1000)
    {
      if (RANK == 0)
        {
          std::cerr << "------- step " << step << " at t = " << absTime << " -------- \n";
        }
      // prediction context
      {
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        noParts = Data.size1();
        noGhosts = GData.size1();
        std::cerr << "Rank " << RANK << ": "
                  << noParts << " particles and " << noGhosts << " ghosts\n";

        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        ghostProxies.resize(noParts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }

        TimeStart = microsec_clock::local_time();
        sphlatch::OctTree BarnesHutTree(theta,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        BarnesHutTree.calcMultipoles();
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

        TimeStart = microsec_clock::local_time();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      }

      dt = 0.000005;
      AccInt.Predictor(dt);
      absTime += dt;
      MemManager.SaveParameter("TIME", absTime, true);

      // correction context
      {
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        noParts = Data.size1();
        noGhosts = GData.size1();
        std::cerr << "Rank " << RANK << ": "
                  << noParts << " particles and " << noGhosts << " ghosts\n";

        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        ghostProxies.resize(noParts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }

        TimeStart = microsec_clock::local_time();
        sphlatch::OctTree BarnesHutTree(theta,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        BarnesHutTree.calcMultipoles();
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

        TimeStart = microsec_clock::local_time();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      }
      AccInt.Corrector(dt);

      if ((step % 10) == 0)
        {
          std::vector<int> outputAttrSet;

          std::string outFilename = "out000000.cdat";
          std::string stepString = boost::lexical_cast<std::string>(step);
          outFilename.replace(outFilename.size() - 5 - stepString.size(),
                              stepString.size(), stepString);

          outputAttrSet += ID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M;
          IOManager.SaveCDAT(outFilename, outputAttrSet);
          if (RANK == 0)
            {
              std::cerr << "saved file " << outFilename << "\n";
            }
        }

      step++;
      Erazer();
    }
  
  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}

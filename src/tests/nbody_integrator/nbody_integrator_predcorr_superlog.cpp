// some defs

#define OOSPH_SINGLE_PRECISION
#define SPHLATCH_SINGLEPREC

//#define OOSPH_LOADBALANCE
#define SPHLATCH_PARALLEL

#include <cstdlib>
#include <iostream>
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
  double logStartTime = MPI_Wtime();

  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  //const size_t SIZE = MPI::COMM_WORLD.Get_size();

  po::options_description Options("Global Options");
  Options.add_options() ("help,h", "Produces this Help blabla...")
  ("input-file,i", po::value<std::string>(), "InputFile")
  ("output-file,o", po::value<std::string>(), "OutputFile");

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

  std::string OutputTag = VMap["output-file"].as<std::string>();
  std::string InputFileName = VMap["input-file"].as<std::string>();
  IOManager.LoadCDAT(InputFileName);
  Erazer();
  
  value_type dt, absTime = 0.; // load from file
  value_type theta = 0.7;
  size_t noParts, noGhosts, step = 0;

  std::string logFilename = "logRank000";
  std::string rankString = boost::lexical_cast<std::string>(RANK);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    start log\n";
  
  // bootstrapping context
  if (RANK == 0)
    {
      std::cerr << "-------- bootstrapping step -----\n";
    }
  {
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << " BB bootstrap\n" << std::flush;
    ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    exchanged parts\n" << std::flush;
    ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    exchanged ghosts\n" << std::flush;
    noParts = Data.size1();
    noGhosts = GData.size1();
    std::cerr << "Rank " << RANK << ": "
              << noParts << " particles and " << noGhosts << " ghosts\n";

    partProxies.resize(noParts);
    for (size_t i = 0; i < noParts; i++)
      {
        (partProxies[i]).setup(&Data, i);
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    prepared part proxies\n" << std::flush;
    ghostProxies.resize(noGhosts);
    for (size_t i = 0; i < noGhosts; i++)
      {
        (ghostProxies[i]).setup(&GData, i);
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    prepared ghost proxies\n" << std::flush;

    TimeStart = microsec_clock::local_time();
    sphlatch::OctTree BarnesHutTree(theta,
                                    1.0,
                                    CostZone.getDepth(),
                                    CostZone.getCenter(),
                                    CostZone.getSidelength());
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    constructed tree\n" << std::flush;
    for (size_t i = 0; i < noParts; i++)
      {
        BarnesHutTree.insertParticle(*(partProxies[i]), true);
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    inserted parts\n" << std::flush;
    for (size_t i = 0; i < noGhosts; i++)
      {
        BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    inserted ghosts\n" << std::flush;

    MPI::COMM_WORLD.Barrier();
    
BarnesHutTree.calcMultipoles();
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    calculated multipoles\n" << std::flush;
    TimeStop = microsec_clock::local_time();
    std::cerr << "Rank " << RANK << ": "
              << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

    TimeStart = microsec_clock::local_time();
    for (size_t i = 0; i < noParts; i++)
      {
        BarnesHutTree.calcGravity(*(partProxies[i]));
      }
    logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    calculated gravity\n" << std::flush;
    TimeStop = microsec_clock::local_time();
    std::cerr << "Rank " << RANK << ": "
              << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
  }
  AccInt.BootStrap();
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    integrator bootstrapped\n" << std::flush;
  Erazer();
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    erazed\n" << std::flush;
  

  while (step < 100000)
    {
      if (RANK == 0)
        {
          std::cerr << "------- step " << step << " at t = " << absTime << " -------- \n";
        }
      // prediction context
      {
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << " PP predicting step " << step << "\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    exchanged parts\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    exchanged ghosts\n" << std::flush;
        noParts = Data.size1();
        noGhosts = GData.size1();
        std::cerr << "Rank " << RANK << ": "
                  << noParts << " particles and " << noGhosts << " ghosts\n";

        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    prepared part proxies\n" << std::flush;
        ghostProxies.resize(noGhosts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    prepared ghost proxies\n" << std::flush;

        TimeStart = microsec_clock::local_time();
        sphlatch::OctTree BarnesHutTree(theta,
                                        1.0,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    constructed tree\n" << std::flush;
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    inserted parts\n" << std::flush;
        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    inserted ghosts\n" << std::flush;
        BarnesHutTree.calcMultipoles();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    calculated multipoles\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

        TimeStart = microsec_clock::local_time();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    calculated gravity\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      }

      //dt = 0.000005;
      dt = 0.000001;
      AccInt.Predictor(dt);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    integrator predicted\n" << std::flush;
      absTime += dt;
      MemManager.SaveParameter("TIME", absTime, true);

      // correction context
      {
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << " CC correcting step " << step << "\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    exchanged parts\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    exchanged ghosts\n" << std::flush;
        noParts = Data.size1();
        noGhosts = GData.size1();
        std::cerr << "Rank " << RANK << ": "
                  << noParts << " particles and " << noGhosts << " ghosts\n";

        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    prepared part proxies\n" << std::flush;
        ghostProxies.resize(noGhosts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    prepared ghost proxies\n" << std::flush;

        TimeStart = microsec_clock::local_time();
        sphlatch::OctTree BarnesHutTree(theta,
                                        1.0,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    constructed tree\n" << std::flush;
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    inserted parts\n" << std::flush;
        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    inserted ghosts\n" << std::flush;
        BarnesHutTree.calcMultipoles();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    calculated multipoles\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

        TimeStart = microsec_clock::local_time();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    calculated gravity\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      }
      AccInt.Corrector(dt);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    integrator corrected\n" << std::flush;
      
      if ((step % 100) == 0)
        {
          std::vector<int> outputAttrSet;

          std::string outFilename = OutputTag;
	  outFilename += "000000.cdat";
          std::string stepString = boost::lexical_cast<std::string>(step);
          outFilename.replace(outFilename.size() - 5 - stepString.size(),
                              stepString.size(), stepString);

          outputAttrSet += ID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M;
          IOManager.SaveCDAT(outFilename, outputAttrSet);
          if (RANK == 0)
            {
              std::cerr << "saved file " << outFilename << "\n";
            }
          logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    saved file\n" << std::flush;
        }

      step++;
      Erazer();
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6) << MPI_Wtime() - logStartTime << "    erazed\n" << std::flush;
    }
  
  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}

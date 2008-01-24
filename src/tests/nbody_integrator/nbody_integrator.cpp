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
typedef SimTrait::index_vector_type index_vector_type;

#include "iomanager.h"
typedef oosph::IOManager<SimTrait> io_type;

#include "memorymanager.h"
typedef oosph::MemoryManager<SimTrait> mem_type;

#include "communicationmanager.h"
typedef oosph::CommunicationManager<SimTrait> com_type;

#include "costzone.h"
typedef mpl::vector_c<size_t, oosph::X> CostZoneIndex;
typedef oosph::CostZone<CostZoneIndex, SimTrait> CostZoneType;

#include "verlet.h"
typedef mpl::vector_c<size_t, oosph::X, oosph::VX, oosph::AX, oosph::OAX> AccIntIndices;
typedef oosph::VecVerlet<AccIntIndices, SimTrait> AccIntType;

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

  po::options_description Options("Global Options");
  Options.add_options() ("help,h", "Produces this Help")
  ("input-file,i", po::value<std::string>(), "input file")
  ("output-tag,o", po::value<std::string>(), "tag for output files");

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

  Data.resize(Data.size1(), oosph::PX); // Verlet only needs vars up to OAZ
  GData.resize(GData.size1(), oosph::OX);

  std::string InputFileName = VMap["input-file"].as<std::string>();
  std::string outputTag = VMap["output-tag"].as<std::string>();
  IOManager.LoadCDAT(InputFileName);
  value_type dt, absTime = 0.; // load from file
  value_type theta = 0.7;
  size_t noParts, noGhosts, step = 0;

  std::string logFilename = "logRank000";
  std::string rankString = boost::lexical_cast<std::string>(RANK);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "    start log\n";

  while (step < 10000)
    {
      if (RANK == 0)
        {
          std::cerr << "------- step " << step << " at t = " << absTime << " -------- \n";
        }

      const value_type maxRadius = 1.;
      index_vector_type delParts;
      for (size_t j = 0; j < Data.size1(); j++)
        {
          value_type curRadius = sqrt(
            (Data(j, X) - 0.5) * (Data(j, X) - 0.5) +
            (Data(j, Y) - 0.5) * (Data(j, Y) - 0.5) +
            (Data(j, Z) - 0.5) * (Data(j, Z) - 0.5));
          if (curRadius > maxRadius)
            {
              delParts += j;
            }
        }
      MemManager.DelParticles(delParts);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime
              << "    deleted " << delParts.size() << " particles\n" << std::flush;

      // next step position context
      {
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime
                << "    predicting step " << step << "\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        noParts = Data.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime
                << "    exchanged parts  (" << noParts << ")\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        noGhosts = GData.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime
                << "    exchanged ghosts (" << noGhosts << ")\n" << std::flush;

        std::cerr << "Rank " << RANK << ": "
                  << noParts << " particles and " << noGhosts << " ghosts\n";

        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    prepared part proxies\n" << std::flush;

        ghostProxies.resize(noGhosts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    prepared ghost proxies\n" << std::flush;

        TimeStart = microsec_clock::local_time();
        sphlatch::OctTree BarnesHutTree(theta,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    constructed tree\n" << std::flush;

        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    inserted parts\n" << std::flush;

        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    inserted ghosts\n" << std::flush;
        
        MPI::COMM_WORLD.Barrier();
        BarnesHutTree.calcMultipoles();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    calculated multipoles\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

        TimeStart = microsec_clock::local_time();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    calculated gravity\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      }

      dt = 0.000001; // just right for random_1M.cdat
      AccInt.getPos(dt);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime
              << "    intergrator next position\n" << std::flush;
      absTime += dt;
      MemManager.SaveParameter("TIME", absTime, true);

      // next step velocity context
      {
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime
                << "    predicting step " << step << "\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        noParts = Data.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime
                << "    exchanged parts  (" << noParts << ")\n" << std::flush;
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        noGhosts = GData.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime
                << "    exchanged ghosts (" << noGhosts << ")\n" << std::flush;

        std::cerr << "Rank " << RANK << ": "
                  << noParts << " particles and " << noGhosts << " ghosts\n";

        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    prepared part proxies\n" << std::flush;

        ghostProxies.resize(noGhosts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    prepared ghost proxies\n" << std::flush;

        TimeStart = microsec_clock::local_time();
        sphlatch::OctTree BarnesHutTree(theta,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    constructed tree\n" << std::flush;

        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    inserted parts\n" << std::flush;

        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    inserted ghosts\n" << std::flush;
                
        MPI::COMM_WORLD.Barrier();
        BarnesHutTree.calcMultipoles();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    calculated multipoles\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "B&H tree build time     " << (TimeStop - TimeStart) << "\n";

        TimeStart = microsec_clock::local_time();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - logStartTime << "    calculated gravity\n" << std::flush;
        TimeStop = microsec_clock::local_time();
        std::cerr << "Rank " << RANK << ": "
                  << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      }
      
      AccInt.getVel(dt);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime
              << "    intergrator next velocities\n" << std::flush;

      if ((step % 10) == 0)
        {
          std::vector<int> outputAttrSet;

          std::string outFilename = outputTag;
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
          logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                  << MPI_Wtime() - logStartTime << "    saved file\n" << std::flush;
        }

      step++;
      Erazer();
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - logStartTime << "    erazed\n" << std::flush;
    }

  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}

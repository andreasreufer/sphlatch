// some defs

// uncomment for single-precision calculation
//#define OOSPH_SINGLE_PRECISION
//#define SPHLATCH_SINGLEPREC

// enable load-balancing (fishy!)
//#define OOSPH_LOADBALANCE

// enable parallel tree
#define SPHLATCH_MPI

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

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
  ("output-tag,o", po::value<std::string>(), "tag for output files")
  ("time-step,d", po::value<value_type>(), "integration timestep dt")
  ("save-time,t", po::value<value_type>(), "save a dump every dt_save")
  ("stop-time,s", po::value<value_type>(), "stop the simulation at stoptime");


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
      !VMap.count("output-tag") or
      !VMap.count("time-step") or
      !VMap.count("save-time") or
      !VMap.count("stop-time"))
    {
      std::cerr << Options << std::endl;
      MPI::Finalize();
      return EXIT_FAILURE;
    }

  io_type& IOManager(io_type::Instance());
  mem_type& MemManager(mem_type::Instance());
  com_type& ComManager(com_type::Instance());
  CostZoneType& CostZone(CostZoneType::Instance());

  ErazerType& Erazer(ErazerType::Instance());
  AccIntType AccInt;

  std::vector<sphlatch::NodeProxy> partProxies, ghostProxies;

  SimTrait::matrix_reference Data(MemManager.Data);
  SimTrait::matrix_reference GData(MemManager.GData);

  Data.resize(Data.size1(), oosph::PX); // Verlet only needs vars up to OAZ
  GData.resize(GData.size1(), oosph::OX);
  size_t noParts, noGhosts;

  std::string inputFileName = VMap["input-file"].as<std::string>();
  std::string outputTag = VMap["output-tag"].as<std::string>();

  IOManager.LoadCDAT(inputFileName);

  value_type gravTheta = MemManager.LoadParameter("GRAVTHETA");
  if (gravTheta != gravTheta)
    {
      gravTheta = 0.7;       // standard value for opening angle theta
    }
  MemManager.SaveParameter("GRAVTHETA", gravTheta, true);

  value_type gravConst = MemManager.LoadParameter("GRAVCONST");
  if (gravConst != gravConst)
    {
      //gravConst = 1.;       //  G=1 system
      gravConst = 6.67259e-11;       //  SI units
    }
  MemManager.SaveParameter("GRAVCONST", gravConst, true);

  value_type absTime = MemManager.LoadParameter("TIME");
  if (absTime != absTime)
    {
      absTime = 0.;       // start at t = 0
    }
  MemManager.SaveParameter("TIME", absTime, true);

  value_type maxRadius = MemManager.LoadParameter("MAXRAD");
  if (maxRadius != maxRadius)
    {
      //maxRadius = 1;
      maxRadius = 6.0e11;  // ~ 4 AU in SI-units
    }
  MemManager.SaveParameter("MAXRAD", maxRadius, true);

  identType starID = static_cast<identType>(MemManager.LoadParameter("STARID"));
  if (starID != starID)
    {
      starID = 1;
    }
  MemManager.SaveParameter("starID", starID, true);

  // set up logging stuff
  std::string logFilename = "logRank000";
  std::string rankString = boost::lexical_cast<std::string>(RANK);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << "    start log\n";

  value_type saveTimestep, stopTime, dtInput, dtSave, dt, lastSavetime;

  saveTimestep = VMap["save-time"].as<value_type>();
  stopTime = VMap["stop-time"].as<value_type>();
  dtInput = VMap["time-step"].as<value_type>();
  lastSavetime = absTime - saveTimestep;

  size_t step = 0;

  while (absTime < stopTime)
    {
      if (RANK == 0)
        {
          std::cout << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                    << MPI_Wtime() - logStartTime
                    << "    step " << step << " at t = "
                    << std::scientific << absTime << "\n";
        }

      double stepStartTime = MPI_Wtime();

      index_vector_type delParts;
      for (size_t j = 0; j < Data.size1(); j++)
        {
          value_type curRadius = sqrt(
            (Data(j, X) - 0.0) * (Data(j, X) - 0.0) +
            (Data(j, Y) - 0.0) * (Data(j, Y) - 0.0) +
            (Data(j, Z) - 0.0) * (Data(j, Z) - 0.0));
          if (curRadius > maxRadius)
            {
              delParts += j;
            }
        }
      MemManager.DelParticles(delParts);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "    deleted " << delParts.size() << " particles\n" << std::flush;

      // next step position context
      {
        stepStartTime = MPI_Wtime();
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        noParts = Data.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime
                << "    exchanged parts  (" << noParts << ")\n" << std::flush;

        stepStartTime = MPI_Wtime();
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        noGhosts = GData.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime
                << "    exchanged ghosts (" << noGhosts << ")\n" << std::flush;

        for (size_t i = 0; i < noParts; i++)
          {
            if (starID == static_cast<identType>(Data(i, ID)))
              {
                Data(i, X) = 0.;
                Data(i, Y) = 0.;
                Data(i, Z) = 0.;
                Data(i, VX) = 0.;
                Data(i, VY) = 0.;
                Data(i, VZ) = 0.;
              }
          }

        stepStartTime = MPI_Wtime();
        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    prepared part proxies\n" << std::flush;

        stepStartTime = MPI_Wtime();
        ghostProxies.resize(noGhosts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    prepared ghost proxies\n" << std::flush;

        stepStartTime = MPI_Wtime();
        sphlatch::OctTree BarnesHutTree(gravTheta,
                                        gravConst,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    constructed tree\n" << std::flush;

        stepStartTime = MPI_Wtime();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    inserted parts\n" << std::flush;

        stepStartTime = MPI_Wtime();
        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    inserted ghosts\n" << std::flush;

        stepStartTime = MPI_Wtime();
        BarnesHutTree.calcMultipoles();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    calculated multipoles\n" << std::flush;

        Erazer();
        stepStartTime = MPI_Wtime();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    calculated gravity\n" << std::flush;
      }

      // determine timestep to next save
      stepStartTime = MPI_Wtime();
      dtSave = (lrint(lastSavetime / saveTimestep)
                + 1) * saveTimestep - absTime;
      dt = std::min(dtInput, dtSave);
      AccInt.getPos(dt);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "    intergrate to next position\n" << std::flush;
      absTime += dt;
      MemManager.SaveParameter("TIME", absTime, true);

      // next step velocity context
      {
        stepStartTime = MPI_Wtime();
        ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
        noParts = Data.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime
                << "    exchanged parts  (" << noParts << ")\n" << std::flush;

        stepStartTime = MPI_Wtime();
        ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
        noGhosts = GData.size1();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime
                << "    exchanged ghosts (" << noGhosts << ")\n" << std::flush;

        for (size_t i = 0; i < noParts; i++)
          {
            if (starID == static_cast<identType>(Data(i, ID)))
              {
                Data(i, X) = 0.;
                Data(i, Y) = 0.;
                Data(i, Z) = 0.;
                Data(i, VX) = 0.;
                Data(i, VY) = 0.;
                Data(i, VZ) = 0.;
              }
          }

        stepStartTime = MPI_Wtime();
        partProxies.resize(noParts);
        for (size_t i = 0; i < noParts; i++)
          {
            (partProxies[i]).setup(&Data, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    prepared part proxies\n" << std::flush;

        stepStartTime = MPI_Wtime();
        ghostProxies.resize(noGhosts);
        for (size_t i = 0; i < noGhosts; i++)
          {
            (ghostProxies[i]).setup(&GData, i);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    prepared ghost proxies\n" << std::flush;

        stepStartTime = MPI_Wtime();
        sphlatch::OctTree BarnesHutTree(gravTheta,
                                        gravConst,
                                        CostZone.getDepth(),
                                        CostZone.getCenter(),
                                        CostZone.getSidelength());
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    constructed tree\n" << std::flush;

        stepStartTime = MPI_Wtime();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.insertParticle(*(partProxies[i]), true);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    inserted parts\n" << std::flush;

        stepStartTime = MPI_Wtime();
        for (size_t i = 0; i < noGhosts; i++)
          {
            BarnesHutTree.insertParticle(*(ghostProxies[i]), false);
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    inserted ghosts\n" << std::flush;

        stepStartTime = MPI_Wtime();
        BarnesHutTree.calcMultipoles();
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    calculated multipoles\n" << std::flush;

        Erazer();
        stepStartTime = MPI_Wtime();
        for (size_t i = 0; i < noParts; i++)
          {
            BarnesHutTree.calcGravity(*(partProxies[i]));
          }
        logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                << MPI_Wtime() - stepStartTime << "    calculated gravity\n" << std::flush;
      }

      stepStartTime = MPI_Wtime();
      AccInt.getVel(dt);
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "    intergrate to next velocities\n" << std::flush;

      if (fabs((lrint(lastSavetime / saveTimestep)
                + 1) * saveTimestep - absTime) < 1e-5 * saveTimestep)
        {
          std::vector<int> outputAttrSet;

          std::string outFilename = outputTag;
          outFilename += "000000.cdat";
          std::string stepString = boost::lexical_cast<std::string>(step);
          outFilename.replace(outFilename.size() - 5 - stepString.size(),
                              stepString.size(), stepString);

          if (RANK == 0)
            {
              std::cout << std::fixed << std::right
                        << std::setw(15) << std::setprecision(6)
                        << MPI_Wtime() - logStartTime
                        << "     save dump to " << outFilename << "\n";
            }

          outputAttrSet += ID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M, GRAVEPS;
          IOManager.SaveCDAT(outFilename, outputAttrSet);

          logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                  << MPI_Wtime() - stepStartTime << "    saved dump\n" << std::flush;
          lastSavetime = absTime;
        }

      step++;
    }

  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}

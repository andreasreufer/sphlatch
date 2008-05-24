// some defs

//#define SPHLATCH_CARTESIAN_XYZ
//#define SPHLATCH_CARTESIAN_YZX
//#define SPHLATCH_CARTESIAN_ZXY
#define SPHLATCH_HILBERT3D

// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

// enable parallel version
#define SPHLATCH_PARALLEL

// enable intensive logging for toptree global summation
//#define SPHLATCH_TREE_LOGSUMUPMP

//#define GRAVITY
#define TREEORDER
//#define RSORDER

//#define SPHLATCH_RANKSPACESERIALIZE

//#define BFCHECK
//#define CHECK_TREE
//#define CHECK_RANKSPACE

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
typedef sphlatch::valueType valueType;

#include "mpi.h"

#include "particle_manager.h"
typedef sphlatch::ParticleManager part_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager com_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

using namespace boost::assign;

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  using namespace boost::assign;
  using namespace sphlatch;

  double stepStartTime, logStartTime = MPI_Wtime();

  com_type& ComManager(com_type::instance());
  part_type& PartManager(part_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  const size_t myDomain = ComManager.getMyDomain();
  const size_t noDomains = ComManager.getNoDomains();

  PartManager.useBasicSPH();
  PartManager.useTimedepH();
  PartManager.useGravity();

  // set up logging stuff
  std::string logFilename = "log_rank000";
  std::string rankString = boost::lexical_cast<std::string>(myDomain);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << " XXX  start log\n";
  
  stepStartTime = MPI_Wtime();

  size_t noParts = 2000000;
  size_t noGhost = 0;
  PartManager.setNoParts(noParts);
  PartManager.resizeAll();

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  for (size_t i = 0; i < noParts; i++)
    {
      id(i) = noParts * myDomain + i + 1;
      m(i) = noParts * myDomain + i + 1;
      h(i) = noParts * myDomain + i + 1;

      noneigh(i) = 50;

      pos(i, X) = static_cast<sphlatch::valueType>(rand()) / RAND_MAX;
      pos(i, Y) = static_cast<sphlatch::valueType>(rand()) / RAND_MAX;
      pos(i, Z) = static_cast<sphlatch::valueType>(rand()) / RAND_MAX;

      vel(i, X) = 0;
      vel(i, Y) = 0;
      vel(i, Z) = 0;
    }

  PartManager.step = 0;
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << "      prepared data (" << noParts << ")\n" << std::flush;

  quantsType exchQuants;
  exchQuants.vects += &pos, &vel;
  exchQuants.ints += &id, &noneigh;
  exchQuants.scalars += &m, &h;

  for (size_t step = 0; step < 65535; step++)
    {
      if ( myDomain == 0 )
      {
        std::cout << std::fixed << std::right << std::setw(15) << std::setprecision(6)
                  << step << "th iteration\n";
      }

      stepStartTime = MPI_Wtime();
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "  X   iteration #" << step << "         "
              << "\n" << std::flush;
      
      CostZone.createDomainPartsIndex();
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      new domain decomposition          "
              << "\n" << std::flush;

      ComManager.exchange(CostZone.domainPartsIndex,
                          CostZone.getNoGhosts(),
                          exchQuants);

      noParts = PartManager.getNoLocalParts();
      noGhost = PartManager.getNoGhostParts();

      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      exchanged data, have now          "
              << noParts << " particles \n" << std::flush;

      ComManager.sendGhostsPrepare(CostZone.createDomainGhostIndex());

      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      prepared ghost exchange expecting "
              << noGhost << " ghosts \n" << std::flush;

      ComManager.sendGhosts(pos);
      ComManager.sendGhosts(vel);
      ComManager.sendGhosts(id);
      ComManager.sendGhosts(m);
      ComManager.sendGhosts(h);

      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      exchanged ghosts \n" << std::flush;

      ///
      /// now move the particle a little bit
      ///
      for (size_t i = 0; i < noParts; i++)
        {
          pos(i, X) += 0.05 * (
            static_cast<sphlatch::valueType>(rand()) / RAND_MAX - 0.5);
          pos(i, Y) += 0.05 * (
            static_cast<sphlatch::valueType>(rand()) / RAND_MAX - 0.5);
          pos(i, Z) += 0.05 * (
            static_cast<sphlatch::valueType>(rand()) / RAND_MAX - 0.5);

          ///
          /// this ensures the particles to stay in a periodic boundary box
          ///
          if ( pos(i, X) > 1.00 )
            pos(i, X) -= 1.00;
          if ( pos(i, Y) > 1.00 )
            pos(i, Y) -= 1.00;
          if ( pos(i, Z) > 1.00 )
            pos(i, Z) -= 1.00;
          
          if ( pos(i, X) < 0.00 )
            pos(i, X) += 1.00;
          if ( pos(i, Y) < 0.00 )
            pos(i, Y) += 1.00;
          if ( pos(i, Z) < 0.00 )
            pos(i, Z) += 1.00;
        }

      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      moved particles by 0.05 at most \n" << std::flush;
    }

  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << " XXX close logfile \n" << std::flush;

  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
};





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

#include "io_manager.h"
typedef sphlatch::IOManager io_type;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager com_type;

#include "costzone.h"
typedef sphlatch::CostZone costzone_type;

using namespace boost::assign;

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

  com_type& ComManager(com_type::instance());
  part_type& PartManager(part_type::instance());
  io_type& IOManager(io_type::instance());
  costzone_type& CostZone(costzone_type::instance());

  const size_t myDomain = ComManager.getMyDomain();
  const size_t noDomains = ComManager.getNoDomains();

  std::string inputFileName = VMap[ "input-file"].as<std::string>();
  std::string outputFileName = VMap["output-file"].as<std::string>();

  PartManager.useBasicSPH();
  PartManager.useTimedepH();
  PartManager.useGravity();

  // set up logging stuff
  std::string logFilename = outputFileName + "_rank000";
  std::string rankString = boost::lexical_cast<std::string>(myDomain);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);
  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - logStartTime << " XXX  start log\n";
  stepStartTime = MPI_Wtime();

  sphlatch::valvectType aux;
  std::cout << "aux var has size " << aux.size() << "\n";
  PartManager.regQuantity(aux, "aux");
  PartManager.regQuantity(aux, "aux1");
  PartManager.setNoParts(500);
  PartManager.resizeAll();
  std::cout << "aux var has size " << aux.size() << "\n";
  PartManager.unRegQuantity(aux);

  /*PartManager.setNoParts(100 + 9*myDomain, 0);
     PartManager.resizeAll();
     //PartManager.setNoParts(02, 80);
     //PartManager.setNoParts(500);
     //PartManager.setNoParts(300000, 500000);

     size_t noParts = PartManager.getNoLocalParts();
     //size_t noGhosts = PartManager.getNoGhostParts();
     logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << "    prepared partManager  (" << noParts << ")\n" << std::flush;

     sphlatch::matrixRefType pos( PartManager.pos );
     sphlatch::matrixRefType vel( PartManager.vel );
     sphlatch::valvectRefType m( PartManager.m );
     sphlatch::valvectRefType h( PartManager.h );
     sphlatch::idvectRefType id( PartManager.id );

     using namespace sphlatch;

     for (size_t i = 0; i < noParts; i++)
     {
     id(i) = 100000*myDomain + i + 1;

     m(i) = 100000*myDomain + i + 1;
     h(i) = 100000*myDomain + i + 1;

     pos(i,X) = i*10 + 0 + 100000*myDomain;
     pos(i,Y) = i*10 + 1 + 100000*myDomain;
     pos(i,Z) = i*10 + 2 + 100000*myDomain;

     vel(i,X) = i*10 + 0 + 100000*myDomain;
     vel(i,Y) = i*10 + 1 + 100000*myDomain;
     vel(i,Z) = i*10 + 2 + 100000*myDomain;
     }

     logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << "      prepare data \n" << std::flush;

     std::set<matrixPtrType> saveVects;
     std::set<valvectPtrType> saveScalars;
     std::set<idvectPtrType> saveInts;

     using namespace boost::assign;
     saveVects += &pos, &vel;
     saveScalars += &m, &h;
     saveInts += &id;

     PartManager.step = 2;
     IOManager.setDoublePrecOut();
     IOManager.setSinglePrecOut();
     IOManager.saveDump( outputFileName, saveVects, saveScalars, saveInts );

     PartManager.attributes["time"] = 3.14;
     PartManager.attributes["gravG"] = 6.0;

     PartManager.step = 29659;
     IOManager.setSinglePrecOut();
     IOManager.saveDump( outputFileName, saveVects, saveScalars, saveInts );

     logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << "      wrote   data \n" << std::flush;

     //std::cout << myDomain << ": " << pos.size1() << "\n";
     //std::cout << pos << "\n";
     //std::cout << h << "\n";

     PartManager.attributes[ "time" ] = 0.;

     PartManager.step = 0;*/

  using namespace sphlatch;

  matrixRefType pos(PartManager.pos);
  matrixRefType vel(PartManager.vel);
  valvectRefType m(PartManager.m);
  valvectRefType h(PartManager.h);
  valvectRefType p(PartManager.p);
  valvectRefType rho(PartManager.rho);
  idvectRefType id(PartManager.id);
  idvectRefType noneigh(PartManager.noneigh);

  stepStartTime = MPI_Wtime();
  IOManager.loadDump(inputFileName);

  //std::cout << PartManager.attributes[ "time" ] << "\n";
  //std::cout << PartManager.step << "\n";
  //std::cout << myDomain << ": " << pos.size1() << "\n";
  //std::cout << pos << "\n";
  //std::cout << h << "\n";

  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << "      read    data \n" << std::flush;

  //std::cout << myDomain << ": " << PartManager.getNoLocalParts() << " parts\n";

  quantsType exchQuants;
  exchQuants.vects += &pos, &vel;
  exchQuants.ints += &id, &noneigh;
  exchQuants.scalars += &m, &h;

  /*for (size_t i = 0; i < 16; i++)
    {
      const size_t noParts = PartManager.getNoLocalParts();

      std::cout << myDomain << ": has " << noParts << " parts\n";

      stepStartTime = MPI_Wtime();

      CostZone.createDomainPartsIndex();
      CostZone.createDomainGhostIndex();

      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      costzone prepared \n" << std::flush;

      stepStartTime = MPI_Wtime();
      ComManager.exchange(CostZone.domainPartsIndex,
                          CostZone.domainGhostIndex,
                          exchQuants);

      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      data exchanged "
              << PartManager.getNoLocalParts() << " "
              << PartManager.getNoGhostParts() << "\n" << std::flush;

      stepStartTime = MPI_Wtime();
      //for (size_t j = 0; j < 10; j++)
      {
      ComManager.sendGhosts( pos );
      ComManager.sendGhosts( vel );
      ComManager.sendGhosts( m  );
      ComManager.sendGhosts( rho);
      ComManager.sendGhosts( id );
      }
      logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
              << MPI_Wtime() - stepStartTime
              << "      ghost pos sent\n" << std::flush;

      for (size_t i = 0; i < noParts; i++)
        {
          pos(i, X) += 0.05 * (
            static_cast<sphlatch::valueType>(rand()) / RAND_MAX - 0.5);
          pos(i, Y) += 0.05 * (
            static_cast<sphlatch::valueType>(rand()) / RAND_MAX - 0.5);
          pos(i, Z) += 0.05 * (
            static_cast<sphlatch::valueType>(rand()) / RAND_MAX - 0.5);
        }
    }*/

  /*
  std::cout << myDomain << ": has " << PartManager.getNoLocalParts() << " particles\n";
  
  CostZone.createDomainPartsIndex();
  
  stepStartTime = MPI_Wtime();
  ComManager.exchange( CostZone.domainPartsIndex,
                       CostZone.getNoGhosts(),
                       exchQuants );

  std::cout << myDomain << ": has " << PartManager.getNoLocalParts() << " particles\n";

  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << "      2nd data xchange\n" << std::flush;

  ComManager.sendGhostsPrepare( CostZone.createDomainGhostIndex() );
  ComManager.sendGhosts(pos);

  std::cout << myDomain << ": sent ghosts!\n";
*/

  domainPartsIndexType domainPartsIndex(2), domainGhostIndex(2);
  
  for (size_t i = 0; i < 2; i++)
  {
    if ( myDomain == i )
    {
      domainPartsIndex[i].resize(20);
      domainGhostIndex[i].resize(0);

      for (size_t j = 0; j < 20; j++)
      {
        (domainPartsIndex[i])[j] = j;
      }
    } else {
      domainPartsIndex[i].resize(10);
      domainGhostIndex[i].resize(10);

      for (size_t j = 0; j < 10; j++)
      {
        (domainPartsIndex[i])[j] = j + 20;
        (domainGhostIndex[i])[j] = j + 10;
      }
    }
  }

  PartManager.setNoParts(30);
  PartManager.resizeAll();
  std::cout << PartManager.getNoLocalParts() << "\n";
  std::cout << pos.size1() << "\n";

  for (size_t i = 0; i < 30; i++)
  {
    pos(i, X) = (myDomain+1)*1000 + i + 1 + 100;
    pos(i, Y) = (myDomain+1)*1000 + i + 1 + 200;
    pos(i, Z) = (myDomain+1)*1000 + i + 1 + 300;
    
    vel(i, X) = (myDomain+1)*1000 + i + 1 + 100;
    vel(i, Y) = (myDomain+1)*1000 + i + 1 + 200;
    vel(i, Z) = (myDomain+1)*1000 + i + 1 + 300;


    m(i) = (myDomain+1);
    h(i) = (myDomain+1);
    
    id(i) = (myDomain+1)*1000 + i;
  }

  if ( myDomain == 0 )
  {
    //std::cout << id << "\n\n";
    std::cout << "no of parts to domain 0    " << domainPartsIndex[0].size() << "\n";
  } else
  {
    std::cout << "no of parts to domain 1    " << domainPartsIndex[0].size() << "\n";
  }

  ComManager.exchange( domainPartsIndex,
                       //domainGhostIndex,
                       10,
                       exchQuants );

  ComManager.sendGhostsPrepare( domainGhostIndex );
  
  /*if ( myDomain == 1 )
  {
    std::cout << "id: \n";
    for (size_t i = 0; i < id.size(); i++)
    {
      std::cout << "  " << i << ":  " << id(i) << "\n";
    }
    std::cout << "\n";
  }*/

  ComManager.sendGhosts( pos );
  ComManager.sendGhosts( id );
  ComManager.sendGhosts( h );

  if ( myDomain == 1 )
  {
    sleep(1);
  }
    std::cout << "id: \n";
    for (size_t i = 0; i < id.size(); i++)
    {
      //std::cout << "  " << i << ":  " << pos(i, X) << "  " << pos(i, Y) << "  " << pos(i, Z) << "\n";
      std::cout << "  " << i << ":  " << id(i) << "\n";
    }
    std::cout << "\n";
  //std::cout << id << "\n";

  logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
          << MPI_Wtime() - stepStartTime
          << " XXX  close logfile \n" << std::flush;

  logFile.close();
  MPI::Finalize();
  return EXIT_SUCCESS;
}



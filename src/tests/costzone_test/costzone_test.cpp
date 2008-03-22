#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
namespace po = boost::program_options;

#include <boost/assign/std/vector.hpp>
using namespace boost::assign;

//#define SPHLATCH_SINGLEPREC
#define SPHLATCH_MPI

#include "particle.h"
#include "iomanager.h"
typedef sphlatch::IOManager io_type;
#include "memorymanager.h"
typedef sphlatch::MemoryManager mem_type;
#include "communicationmanager.h"
typedef sphlatch::CommunicationManager comm_type;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
  MPI::Init(argc, argv);
#endif
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
      return EXIT_FAILURE;
    }

  if (!VMap.count("input-file"))
    {
      std::cout << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type& IOManager(io_type::instance());
  mem_type& MemManager(mem_type::instance());
  comm_type& CommManager(comm_type::instance());

  sphlatch::matrixRefType Data(MemManager.Data);

  std::string InputFileName = VMap["input-file"].as<std::string>();

  Data.resize(Data.size1(), sphlatch::SIZE);
  IOManager.loadCDAT(InputFileName);

  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  const size_t SIZE = MPI::COMM_WORLD.Get_size();
  const size_t noParts = Data.size1();

  sphlatch::domainPartsIndexType partsDomainMapping;

  size_t sendtoRank = (RANK + 1) % SIZE; // send to next rank

  partsDomainMapping.resize(SIZE);

  for (size_t i = 0; i < noParts; i++)
    {
      (partsDomainMapping[sendtoRank]).push_back(i);
    }


  // particles are all distributed now
  using namespace boost::posix_time;
  ptime TimeStart, TimeStop;
  // set up logging stuff
  std::string logFilename = "logRank000";
  std::string rankString = boost::lexical_cast<std::string>(RANK);
  logFilename.replace(logFilename.size() - 0 - rankString.size(),
                      rankString.size(), rankString);

  std::fstream logFile;
  logFile.open(logFilename.c_str(), std::ios::out);
  //logFile << std::fixed << std::right << std::setw(15) << std::setprecision(6)
  //        << MPI_Wtime() - logStartTime << "    start log\n";

  std::cout << "RANK " << RANK << ": " << Data(0,0) << "   " << Data(0,1) << "\n\n\n";
  for (size_t i = 0; i < SIZE; i++)
    {
      TimeStart = microsec_clock::local_time();
      CommManager.exchange(Data, partsDomainMapping, Data);
      TimeStop = microsec_clock::local_time();
      logFile << "communication " << (TimeStop - TimeStart) << "\n";
      std::cout << "RANK " << RANK << ": " << Data(0,0) << "   " << Data(0,1) << "\n";
    }

  sphlatch::countsVectType counts;
  counts.resize(32768);

  long int locCounter = 0;
  for (size_t i = 0; i < counts.size(); i++)
  {
    counts[i] = ( rand() / 32768 ) % ( 7*(RANK+3) ) ;
    locCounter += counts[i];
  }
  std::cout << "RANK " << RANK << ": " << locCounter << "\n";
  
  TimeStart = microsec_clock::local_time();
  CommManager.sumUpCounts(counts);
  TimeStop = microsec_clock::local_time();
  logFile << "summing up " << (TimeStop - TimeStart) << "\n";

  locCounter = 0;
  for (size_t i = 0; i < counts.size(); i++)
  {
    locCounter += counts[i];
  }
  std::cout << "RANK " << RANK << ": " << locCounter << "\n";

  TimeStart = microsec_clock::local_time();
  if ( RANK % 2 == 0 )
  {
    CommManager.sendMatrix(Data, RANK+1);
  } else {
    CommManager.recvMatrix(Data, RANK-1);
  }
  TimeStop = microsec_clock::local_time();
  logFile << "exchange matrices " << (TimeStop - TimeStart) << "\n";
  

  std::vector<long double> intVect(500000);
  TimeStart = microsec_clock::local_time();
  if ( RANK % 2 == 0 )
  {
    CommManager.sendVector<long double>(intVect, RANK+1);
  } else {
    CommManager.recvVector<long double>(intVect, RANK-1);
  }
  TimeStop = microsec_clock::local_time();
  logFile << "exchange vects " << (TimeStop - TimeStart) << "\n";
  
  sphlatch::bitsetType myBitset(3000000);
  TimeStart = microsec_clock::local_time();
  if ( RANK % 2 == 0 )
  {
    CommManager.sendBitset(myBitset, RANK+1);
  } else {
    CommManager.recvBitset(myBitset, RANK-1);
  }
  TimeStop = microsec_clock::local_time();
  logFile << "exchange bitsets " << (TimeStop - TimeStart) << "\n";

  logFile.close();

  sphlatch::valueType myRank = RANK;
  CommManager.sum(myRank);
  std::cout << myRank << "\n";
  
  using namespace sphlatch;
  std::vector<int> outputAttrSet;
  outputAttrSet += ID, X, Y, Z, AX, AY, AZ, M;
  IOManager.saveCDAT("out.cdat", outputAttrSet);

  #ifdef SPHLATCH_MPI
  MPI::Finalize();
  #endif

  return EXIT_SUCCESS;
}

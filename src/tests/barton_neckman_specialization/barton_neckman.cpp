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

#define SPHLATCH_SINGLEPREC
//#define SPHLATCH_MPI

#include "particle.h"
#include "iomanager.h"
typedef sphlatch::IOManager io_type;
#include "memorymanager.h"
typedef sphlatch::MemoryManager mem_type;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>

#include "bhtree.h"

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

  sphlatch::matrixRefType Data(MemManager.Data);

  std::string InputFileName = VMap["input-file"].as<std::string>();
  
  Data.resize(Data.size1(), sphlatch::SIZE);
  IOManager.loadCDAT(InputFileName);
  
  const size_t noParts = Data.size1();

  // particles are all distributed now
  using namespace boost::posix_time;
  ptime TimeStart, TimeStop;

  std::vector<sphlatch::particleProxy> partProxies;
  partProxies.resize(noParts);
  for (size_t i = 0; i < noParts; i++)
    {
      (partProxies[i]).setup(&Data, i);
    }
  sphlatch::valvectType universeCenter(3);
  universeCenter(0) = 0.0;
  universeCenter(1) = 0.0;
  universeCenter(2) = 0.0;
  sphlatch::valueType universeSize = 10., theta = 0.51;
  size_t costzoneDepth = 0;

  //for (size_t i = 0; i < 256; i++)
  for (size_t i = 0; i < 16; i++)
  //for (size_t i = 0; i < 1; i++)
    {
      TimeStart = microsec_clock::local_time();
      //sphlatch::BHtree<sphlatch::Monopoles> BarnesHutTree(theta, 1.0,
      //sphlatch::BHtree<sphlatch::Quadrupoles> BarnesHutTree(theta, 1.0,
      sphlatch::BHtree<sphlatch::Octupoles> BarnesHutTree(theta, 1.0,
                                                          costzoneDepth,
                                                          universeCenter,
                                                          universeSize);

      TimeStop = microsec_clock::local_time();
      std::cerr << "Tree prepare time       " << (TimeStop - TimeStart) << "\n";

      TimeStart = microsec_clock::local_time();
      bool locality = true;
      for (size_t i = 0; i < noParts; i++)
        {
          BarnesHutTree.insertParticle(*(partProxies[i]), locality);
        }
      TimeStop = microsec_clock::local_time();
      std::cerr << "Tree populate time      " << (TimeStop - TimeStart) << "\n";

      TimeStart = microsec_clock::local_time();
      BarnesHutTree.calcMultipoles();
      TimeStop = microsec_clock::local_time();
      std::cerr << "Calc. multipoles time   " << (TimeStop - TimeStart) << "\n";

      //BarnesHutTree.treeDump("dump.txt");
      //BarnesHutTree.treeDOTDump("dump.dot");

      //boost::progress_display show_progress(noParts, std::cout);
      TimeStart = microsec_clock::local_time();
      for (size_t i = 0; i < noParts; i++)
        {
          BarnesHutTree.calcGravity(*(partProxies[i]));
        }
      TimeStop = microsec_clock::local_time();
      std::cerr << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
      std::cerr << "\n";

  }

  using namespace sphlatch;
  std::vector<int> outputAttrSet;
  outputAttrSet += ID, X, Y, Z, AX, AY, AZ, M;
  IOManager.saveCDAT("out.cdat", outputAttrSet);

  #ifdef SPHLATCH_MPI
  MPI::Finalize();
  #endif

  return EXIT_SUCCESS;
}

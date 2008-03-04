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

  std::string InputFileName = "random.cdat";

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
  sphlatch::valueType universeSize = 10., theta = 0.60;
  size_t costzoneDepth = 4;

  TimeStart = microsec_clock::local_time();
  //sphlatch::BHtree<sphlatch::Monopoles> BarnesHutTree(theta, 1.0,
  sphlatch::BHtree<sphlatch::Quadrupoles> BarnesHutTree(theta, 1.0,
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

  boost::progress_display show_progress(noParts, std::cout);
  TimeStart = microsec_clock::local_time();
  for (size_t i = 0; i < noParts; i++)
    {
      BarnesHutTree.calcGravity(*(partProxies[i]));
      ++show_progress;
    }
  TimeStop = microsec_clock::local_time();
  std::cerr << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";
  std::cerr << "\n";

  using namespace sphlatch;
  std::vector<int> outputAttrSet;
  outputAttrSet += ID, X, Y, Z, AX, AY, AZ, M;
  IOManager.saveCDAT("outBH.cdat", outputAttrSet);

  show_progress.restart(noParts * (noParts - 1));
  TimeStart = microsec_clock::local_time();

  for (size_t i = 0; i < noParts; i++)
    {
      sphlatch::valueType partDistPow2, partDistPow3, epsilonPow2;
      Data(i, AX) = 0.;
      Data(i, AY) = 0.;
      Data(i, AZ) = 0.;

      for (size_t j = 0; j < noParts; j++)
        {
          if (i != j)
            {
              partDistPow2 = (Data(i, X) - Data(j, X)) *
                             (Data(i, X) - Data(j, X)) +
                             (Data(i, Y) - Data(j, Y)) *
                             (Data(i, Y) - Data(j, Y)) +
                             (Data(i, Z) - Data(j, Z)) *
                             (Data(i, Z) - Data(j, Z));
              epsilonPow2 = Data(i, GRAVEPS) * Data(i, GRAVEPS);
              partDistPow3 = pow(partDistPow2 + epsilonPow2, 3. / 2.);

              Data(i, AX) -= Data(j, M) * (Data(i, X) - Data(j, X))
                             / partDistPow3;
              Data(i, AY) -= Data(j, M) * (Data(i, Y) - Data(j, Y))
                             / partDistPow3;
              Data(i, AZ) -= Data(j, M) * (Data(i, Z) - Data(j, Z))
                             / partDistPow3;
              ++show_progress;
            }
        }
    }
  TimeStop = microsec_clock::local_time();
  std::cerr << "Gravity BF calc time    " << (TimeStop - TimeStart) << "\n";

  // save particles
  IOManager.saveCDAT("outBF.cdat", outputAttrSet);

  return EXIT_SUCCESS;
}

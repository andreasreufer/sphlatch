#include <cstdlib>
#include <iostream>
#include <string>

#define OOSPH_STANDALONE
#define OOSPH_SINGLE_PRECISION
#define SPHLATCH_SINGLEPREC

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

using namespace oosph;
using namespace boost::assign;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <vector>

// tree stuff

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

  io_type& IOManager(io_type::Instance());
  mem_type& MemManager(mem_type::Instance());

  SimTrait::matrix_reference Data(MemManager.Data);

  std::string InputFileName = VMap["input-file"].as<std::string>();

  Data.resize(Data.size1(), oosph::SIZE);

  IOManager.LoadCDAT(InputFileName);

  const size_t noParts = Data.size1();

  // particles are all distributed now
  using namespace boost::posix_time;
  ptime TimeStart, TimeStop;

  std::vector<sphlatch::NodeProxy> partProxies;

  partProxies.resize(noParts);
  for (size_t i = 0; i < noParts; i++)
    {
      (partProxies[i]).setup(&Data, i);
    }

  TimeStart = microsec_clock::local_time();
  
  valvectType universeCenter(3);
  universeCenter(0) = 0.5;
  universeCenter(1) = 0.5;
  universeCenter(2) = 0.5;
  
  valueType universeSize = 1., theta = 0.6;
  size_t costzoneDepth = 4;
  
  sphlatch::OctTree BarnesHutTree(theta, 1.0, costzoneDepth, universeCenter, universeSize);
  
  TimeStop = microsec_clock::local_time();
  std::cerr << "Tree prepare time       " << (TimeStop - TimeStart) << "\n";

  TimeStart = microsec_clock::local_time();
  for (size_t i = 0; i < noParts; i++)
    {
      BarnesHutTree.insertParticle(*(partProxies[i]), true);
    }
  TimeStop = microsec_clock::local_time();
  std::cerr << "Tree populate time      " << (TimeStop - TimeStart) << "\n";

  std::string treeDumpFilename;

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
  std::cout << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";

  //const value_type epsilonSquare = 0.0025;
  const value_type epsilonSquare = 0.0000;

  show_progress.restart(noParts * (noParts - 1));
  TimeStart = microsec_clock::local_time();

  for (size_t i = 0; i < noParts; i++)
    {
      value_type partDist, partDistPow3;
      Data(i, AX) = 0.;
      Data(i, AY) = 0.;
      Data(i, AZ) = 0.;

      for (size_t j = 0; j < noParts; j++)
        {
          if (i != j)
            {
              partDist = sqrt((Data(i, X) - Data(j, X)) *
                              (Data(i, X) - Data(j, X)) +
                              (Data(i, Y) - Data(j, Y)) *
                              (Data(i, Y) - Data(j, Y)) +
                              (Data(i, Z) - Data(j, Z)) *
                              (Data(i, Z) - Data(j, Z)));
              // partDistPow3 = partDist*partDist*partDist; // unsoftened
              partDistPow3 = partDist * partDist * partDist +
                             partDist * epsilonSquare;

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
  std::vector<int> outputAttrSet;
  std::string outFilename;
  outFilename += "outBF.cdat";
  outputAttrSet += ID, X, Y, Z, AX, AY, AZ, M;
  IOManager.SaveCDAT(outFilename, outputAttrSet);

  return EXIT_SUCCESS;
}

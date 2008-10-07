// uncomment for single-precision calculation
#define SPHLATCH_SINGLEPREC

#include <iostream>
#include <iomanip>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
namespace po = boost::program_options;

#include <boost/assign/std/vector.hpp>

#include "particle.h"

#include "typedefs.h"
//typedef sphlatch::valueType valueType;

#include "iomanager.h"
typedef sphlatch::IOManager io_type;

#include "memorymanager.h"
typedef sphlatch::MemoryManager mem_type;

using namespace boost::assign;

#include <boost/progress.hpp>
#include <vector>
#include <cmath>

// tree stuff
#include "bhtree.h"

#include "ranklist.h"
#include "rankspace.h"

#include "timer.h"

int main(int argc, char* argv[])
{
  po::options_description Options("Global Options");
  Options.add_options() ("help,h", "Produces this Help")
  ("output-file,o", po::value<std::string>(), "output file");

  po::positional_options_description POD;
  POD.add("input-file", 1);

  po::variables_map VMap;
  po::store(po::command_line_parser(argc, argv).options(Options).positional(POD).run(), VMap);
  po::notify(VMap);

  if (VMap.count("help"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  if (!VMap.count("output-file"))
    {
      std::cerr << Options << std::endl;
      return EXIT_FAILURE;
    }

  io_type& IOManager(io_type::instance());
  mem_type& MemManager(mem_type::instance());

  sphlatch::matrixRefType Data(MemManager.Data);
  Data.resize(Data.size1(), sphlatch::SIZE);

  std::string outputFileName = VMap["output-file"].as<std::string>();

  using namespace sphlatch;
  
  size_t noBlobParts = 1000;

  matrixType blob;
  blob.resize(noBlobParts, SIZE);
  
  for (size_t i = 0; i < noBlobParts; i++)
  {
    blob(i, X) = (static_cast<valueType>(random()) / RAND_MAX ) - 0.5;
    blob(i, Y) = (static_cast<valueType>(random()) / RAND_MAX ) - 0.5;
    blob(i, Z) = (static_cast<valueType>(random()) / RAND_MAX ) - 0.5;
  }

  const size_t noBlobs = 1000;

  Data.resize(noBlobParts*noBlobs, SIZE);

  for (size_t i = 0; i < noBlobs; i++)
  {
    
    const valueType phi = 2.*M_PI
      *(static_cast<valueType>(random()) / RAND_MAX );
    const valueType theta = M_PI
      *(static_cast<valueType>(random()) / RAND_MAX );
    const valueType rad = 20.
      *(static_cast<valueType>(random()) / RAND_MAX );
    //const valueType rad = 20*( 1. / ( sin(M_PI - theta) + 0.2 ) );
    //const valueType rad = 20*( 1. / ( (1 - theta/M_PI) + 0.001) );

    /*const valueType xCen = 1.80*rad*sin(theta)*cos(phi);
    const valueType yCen = 1.80*rad*sin(theta)*sin(phi);*/
    //const valueType zCen = 1.00*rad*cos(theta);
    
    /*const valueType gamma = 0.97*M_PI
      *(static_cast<valueType>(random()) / RAND_MAX - 0.5);
    const valueType zCen = 2.0*tan(gamma);*/
    
    /*const valueType alpha = 0.95*M_PI
      *(static_cast<valueType>(random()) / RAND_MAX - 0.5);
    const valueType xCen = 5.*tan(alpha);
    
    const valueType beta = 0.95*M_PI
      *(static_cast<valueType>(random()) / RAND_MAX - 0.5);
    const valueType yCen = 5.*tan(beta);*/

    const valueType xCen = 4.0*rad*sin(theta)*cos(phi);
    const valueType yCen = 4.0*rad*sin(theta)*sin(phi);
    
    valueType zCen;
    /*if ( (static_cast<valueType>(random()) / RAND_MAX ) > 0.97 )
    {
      zCen = -0.6;
    }
    else
    {
      zCen = 0.6;
    }*/
    zCen = 0.0;

    for (size_t j = 0; j < noBlobParts; j++)
    {
      Data(i*noBlobParts + j, ID) = i*noBlobParts + j;
      Data(i*noBlobParts + j, X) = blob(j, X) + xCen;
      Data(i*noBlobParts + j, Y) = blob(j, Y) + yCen;
      Data(i*noBlobParts + j, Z) = blob(j, Z) + zCen;

      Data(i*noBlobParts + j, VX) = 0.;
      Data(i*noBlobParts + j, VY) = 0.;
      Data(i*noBlobParts + j, VZ) = 0.;
      
      Data(i*noBlobParts + j, M) = 1.;
      Data(i*noBlobParts + j, H) = 1.;
    }
  }


  MemManager.simName = "blobs";
  std::vector<int> outputAttrSet;
  outputAttrSet += ID, X, Y, Z, VX, VY, VZ, M, H;
  IOManager.saveCDAT(outputFileName, outputAttrSet);

  return EXIT_SUCCESS;
}



#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#define SPHLATCH_SINGLEPREC
#define SPHLATCH_PARALLEL

#include "bhtree.h"
#include "particle.h"

//#define NPARTS	0
//#define NPARTS	20
#define NPARTS  50
//#define NPARTS	1000
//#define NPARTS	10000
//#define NPARTS	100000
//#define NPARTS	300000

int main(int argc, char* argv[])
{
  MPI::Init(argc, argv);

  const size_t SIZE = MPI::COMM_WORLD.Get_size();
  const size_t RANK = MPI::COMM_WORLD.Get_rank();
  
  using namespace sphlatch;
  using namespace boost::posix_time;
  
  ptime TimeStart, TimeStop;

  matrixType Data(NPARTS, SIZE);

  valueType CenterX = 0;
  valueType CenterY = 0;
  for (size_t i = 0; i < NPARTS; i++)
    {
      Data(i, ID) = i;
      Data(i, X) = ((valueType)rand()) / RAND_MAX;
      Data(i, Y) = ((valueType)rand()) / RAND_MAX;
      //Data(i, Z) = ( (valueType)rand() ) / RAND_MAX;
      Data(i, Z) = 0.;
      Data(i, M) = 1.;
      CenterX += Data(i, X);
      CenterY += Data(i, Y);
    }

  std::vector<sphlatch::particleProxy> partProxies;
  partProxies.resize(NPARTS);
  for (size_t i = 0; i < NPARTS; i++)
    {
      (partProxies[i]).setup(&Data, i);
    }
  sphlatch::valvectType universeCenter(3);
  universeCenter(0) = 0.5;
  universeCenter(1) = 0.5;
  universeCenter(2) = 0.5;
  sphlatch::valueType universeSize = 1., theta = 0.60;
  size_t costzoneDepth = 4;

  // start tree context
  TimeStart = microsec_clock::local_time();
  sphlatch::BHtree<sphlatch::Monopoles> BarnesHutTree(theta, 1.0,
  //sphlatch::BHtree<sphlatch::Quadrupoles> BarnesHutTree(theta, 1.0,
                                                        costzoneDepth,
                                                        universeCenter,
                                                        universeSize);

  TimeStop = microsec_clock::local_time();
  std::cout << "Tree prepare time       " << (TimeStop - TimeStart) << "\n";

  //boost::progress_display show_progress( NPARTS , std::cout);
  TimeStart = microsec_clock::local_time();
  for (size_t i = 0; i < NPARTS; i++)
    {
      /*if ( RANK == 0 ) {
          if ( Data(i, sphlatch::Z) < 0.50 ) {
              BarnesHutTree.insertParticle( *(DataProxies[i]),true );
          }
          if ( Data(i, sphlatch::Z) > 0.50 && Data(i, sphlatch::Z) < 0.75 ) {
              BarnesHutTree.insertParticle( *(DataProxies[i]),false );
          }
         } else {
          if ( Data(i, sphlatch::Z) > 0.50 ) {
              BarnesHutTree.insertParticle( *(DataProxies[i]),true );
          }
          if ( Data(i, sphlatch::Z) < 0.50 && Data(i, sphlatch::Z) > 0.25 ) {
              BarnesHutTree.insertParticle( *(DataProxies[i]),false );
          }
         }*/
      //BarnesHutTree.insertParticle( *(DataProxies[i]), locality);
      BarnesHutTree.insertParticle(*(partProxies[i]), true);
      //++show_progress;
    }
  TimeStop = microsec_clock::local_time();
  std::cout << "Tree populate time      " << (TimeStop - TimeStart) << "\n";

  // dump tree
  std::string treeDumpFilename;
  treeDumpFilename = "beforeMP_";
  treeDumpFilename += boost::lexical_cast<std::string>(RANK + 1);
  treeDumpFilename += "_of_";
  treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
  treeDumpFilename += ".dot";
  BarnesHutTree.treeDOTDump(treeDumpFilename);

  TimeStart = microsec_clock::local_time();
  BarnesHutTree.calcMultipoles();
  TimeStop = microsec_clock::local_time();
  std::cerr << "Calc. multipoles time   " << (TimeStop - TimeStart) << "\n";


  treeDumpFilename = "afterMP_";
  treeDumpFilename += boost::lexical_cast<std::string>(RANK + 1);
  treeDumpFilename += "_of_";
  treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
  treeDumpFilename += ".dot";
  BarnesHutTree.treeDOTDump(treeDumpFilename);

  // dump toptree
  /*std::string toptreeDumpFilename = "toptree_postsync_";
  toptreeDumpFilename += boost::lexical_cast<std::string>(RANK + 1);
  toptreeDumpFilename += "_of_";
  toptreeDumpFilename += boost::lexical_cast<std::string>(SIZE);
  toptreeDumpFilename += ".txt";
  BarnesHutTree.toptreeDump(toptreeDumpFilename);*/

  //show_progress.restart(NPARTS);
  TimeStart = microsec_clock::local_time();
  for (size_t i = 0; i < NPARTS; i++)
    {
      if (((RANK == 0) && Data(i, sphlatch::Z) < 0.50) ||
          ((RANK == 1) && Data(i, sphlatch::Z) > 0.50))
        {
          BarnesHutTree.calcGravity(*(partProxies[i]));
        }
      //++show_progress;
    }

  // dump tree
  treeDumpFilename = "treedump_";
  treeDumpFilename += boost::lexical_cast<std::string>(RANK + 1);
  treeDumpFilename += "_of_";
  treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
  treeDumpFilename += ".txt";
  BarnesHutTree.treeDump(treeDumpFilename);

  TimeStop = microsec_clock::local_time();
  std::cout << "Gravity calc time       " << (TimeStop - TimeStart) << "\n";

  MPI::Finalize();
  return EXIT_SUCCESS;
}

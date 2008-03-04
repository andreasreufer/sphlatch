#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#define SPHLATCH_SINGLEPREC

#include "bhtree.h"

//#define NPARTS	20
#define NPARTS	100000

int main() {
	//typedef OctTree<size_t> Tree;
		
	using namespace boost::posix_time;
    using namespace sphlatch;
    
	ptime TimeStart, TimeStop;
	
	matrixType Data(NPARTS, SIZE);

	valueType CenterX = 0;
	valueType CenterY = 0;
	for (size_t i = 0; i < NPARTS; i++) {
		Data(i, PID) = i;
		Data(i, X) = ( (valueType)rand() ) / RAND_MAX;
		Data(i, Y) = ( (valueType)rand() ) / RAND_MAX;
		Data(i, Z) = ( (valueType)rand() ) / RAND_MAX;
		//Data(i, Z) = 0.;
		Data(i, M) = 1.;
		CenterX += Data(i, X);
		CenterY += Data(i, Y);
	}
		
	std::vector<sphlatch::particleProxy> partProxies;
	partProxies.resize(NPARTS);
	
	for (size_t i = 0; i < NPARTS; i++) {
		( partProxies[i] ).setup(&Data, i);
	}

	// start tree context
	{  
  valvectType universeCenter(3);
  universeCenter(0) = 0.5;
  universeCenter(1) = 0.5;
  universeCenter(2) = 0.5;
  sphlatch::valueType universeSize = 1., theta = 0.60;
  size_t costzoneDepth = 4;

  TimeStart = microsec_clock::local_time();
  sphlatch::BHtree<sphlatch::Monopoles> BarnesHutTree(theta, 1.0,
  //sphlatch::BHtree<sphlatch::Quadrupoles> BarnesHutTree(theta, 1.0,
                                                        costzoneDepth,
                                                        universeCenter,
                                                        universeSize);

  TimeStop = microsec_clock::local_time();
  std::cout << "Tree prepare time       " << (TimeStop - TimeStart) << "\n";
	
	boost::progress_display show_progress( NPARTS , std::cerr);
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < NPARTS; i++) {
		BarnesHutTree.insertParticle( *(partProxies[i]), true);
		++show_progress;
	}
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Tree populate time      " << ( TimeStop - TimeStart ) << "\n";
	
	TimeStart = microsec_clock::local_time();	
	BarnesHutTree.calcMultipoles();
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Calc. multipoles time   " << ( TimeStop - TimeStart ) << "\n";
	
	show_progress.restart(NPARTS);
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < NPARTS; i++) {
		BarnesHutTree.calcGravity(*(partProxies[i]) );
		++show_progress;
	}
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Gravity calc time       " << ( TimeStop - TimeStart ) << "\n";
	
	TimeStart = microsec_clock::local_time();
	}	// tree gets deleted here!
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Tree delete time        " << ( TimeStop - TimeStart ) << "\n";

	sleep(1);
	
	return EXIT_SUCCESS;
}

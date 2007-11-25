#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#define SINGLEPREC

#include "octree.h"
#include "nodeproxy.h"
typedef NodeProxy	NodeProxyType;
typedef NodeProxy*	NodeProxyPtrType;
typedef NodeProxy&	NodeProxyRefType;

//#define NPARTS	20
#define NPARTS	100000

int main() {
	//typedef OctTree<size_t> Tree;
		
	using namespace boost::posix_time;
	ptime TimeStart, TimeStop;
	
	matrixType Data(NPARTS, PSIZE);

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

	/*std::cout << "CenterX " << CenterX / NPARTS
			  << "       CenterY " << CenterY / NPARTS << "\n";*/
		
	std::vector<NodeProxy> DataProxies;
	DataProxies.resize(NPARTS);
	
	for (size_t i = 0; i < NPARTS; i++) {
		( DataProxies[i] ).setup(&Data, i);
	}

	// start tree context
	{
	TimeStart = microsec_clock::local_time();	
//	OctTree<NodeProxyPtrType> BarnesHutTree;
	OctTree BarnesHutTree;

	TimeStop  = microsec_clock::local_time();
	std::cerr << "Tree prepare time       " << ( TimeStop - TimeStart ) << "\n";
	
	boost::progress_display show_progress( NPARTS , std::cerr);
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < NPARTS; i++) {
		BarnesHutTree.insertParticle( *(DataProxies[i]), true);
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
		BarnesHutTree.calcGravity(*(DataProxies[i]) );
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

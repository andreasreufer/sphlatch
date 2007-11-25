#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "octree.h"
#include "nodeproxy.h"
typedef NodeProxy	NodeProxyType;
typedef NodeProxy*	NodeProxyPtrType;
typedef NodeProxy&	NodeProxyRefType;

//#define NPARTS	20
#define NPARTS	30000

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
	TimeStart = microsec_clock::local_time();	
//	OctTree<NodeProxyPtrType> BarnesHutTree;
	OctTree BarnesHutTree;

	TimeStop  = microsec_clock::local_time();
	std::cerr << "Tree prepare time       " << ( TimeStop - TimeStart ) << "\n";
	
	std::vector<NodeProxy> DataProxies;
	DataProxies.resize(NPARTS);
	
	for (size_t i = 0; i < NPARTS; i++) {
		( DataProxies[i] ).setup(&Data, i);
	}
	
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
	
	matrixType TreeData = Data;
	
	show_progress.restart(NPARTS*(NPARTS-1));
	
	TimeStart = microsec_clock::local_time();
	for (size_t i = 0; i < NPARTS; i++) {
		valueType partDist, partDistPow3;
		Data(i, AX) = 0.;
		Data(i, AY) = 0.;
		Data(i, AZ) = 0.;
		
		for (size_t j = 0; j < NPARTS; j++) {
			if ( i != j ) {
				partDist = sqrt(	( Data(i, X) - Data(j, X) )*
									( Data(i, X) - Data(j, X) ) +
									( Data(i, Y) - Data(j, Y) )*
									( Data(i, Y) - Data(j, Y) ) +
									( Data(i, Z) - Data(j, Z) )*
									( Data(i, Z) - Data(j, Z) ) );
				partDistPow3 = partDist*partDist*partDist;
				Data(i, AX) -= Data(j, M) * ( Data(i, X) - Data(j, X) )
								/ partDistPow3;
				Data(i, AY) -= Data(j, M) * ( Data(i, Y) - Data(j, Y) )
								/ partDistPow3;
				Data(i, AZ) -= Data(j, M) * ( Data(i, Z) - Data(j, Z) )
								/ partDistPow3;
				++show_progress;
			}
		}
	}
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Gravity BF calc time    " << ( TimeStop - TimeStart ) << "\n";

	matrixType BFData = Data;
	
	std::vector<valueType> accLengthRelErr;
	accLengthRelErr.resize(NPARTS);
	
	for (size_t i = 0; i < NPARTS; i++) {
		valueType bfAccLength, treeAccLength;
		
		bfAccLength = sqrt(	BFData(i, AX) * BFData(i, AX)  +
							BFData(i, AY) * BFData(i, AY)  +
							BFData(i, AZ) * BFData(i, AZ) );

		treeAccLength = sqrt(	TreeData(i, AX) * TreeData(i, AX) +
								TreeData(i, AY) * TreeData(i, AY) +
								TreeData(i, AZ) * TreeData(i, AZ) );
							
		accLengthRelErr[i] = ( treeAccLength - bfAccLength ) / bfAccLength;
		std::cout << accLengthRelErr[i] << "\n";
	}
	std::cout << "\n";
	
	return EXIT_SUCCESS;
}

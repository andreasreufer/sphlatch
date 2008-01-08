#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//#define SINGLEPREC
#define OOSPH_MPI
#include <mpi.h>

#include "typedefs.h"
#include "octree.h"
#include "nodeproxy.h"
/*typedef NodeProxy	NodeProxyType;
typedef NodeProxy*	NodeProxyPtrType;
typedef NodeProxy&	NodeProxyRefType;*/

//#define NPARTS	20
#define NPARTS	10000

int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);

	//typedef OctTree<size_t> Tree;
		
	using namespace boost::posix_time;
	ptime TimeStart, TimeStop;
	
	matrixType Data(NPARTS, knack::PSIZE);

	valueType CenterX = 0;
	valueType CenterY = 0;
	for (size_t i = 0; i < NPARTS; i++) {
		Data(i, knack::PID) = i;
		Data(i, knack::X) = ( (valueType)rand() ) / RAND_MAX;
		Data(i, knack::Y) = ( (valueType)rand() ) / RAND_MAX;
		Data(i, knack::Z) = ( (valueType)rand() ) / RAND_MAX;
		//Data(i, Z) = 0.;
		Data(i, knack::M) = 1.;
		CenterX += Data(i, knack::X);
		CenterY += Data(i, knack::Y);
	}

	/*std::cout << "CenterX " << CenterX / NPARTS
			  << "       CenterY " << CenterY / NPARTS << "\n";*/
	TimeStart = microsec_clock::local_time();	
//	OctTree<NodeProxyPtrType> BarnesHutTree;
	knack::OctTree BarnesHutTree;

	TimeStop  = microsec_clock::local_time();
	std::cerr << "Tree prepare time       " << ( TimeStop - TimeStart ) << "\n";
	
	std::vector<knack::NodeProxy> DataProxies;
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
		Data(i, knack::AX) = 0.;
		Data(i, knack::AY) = 0.;
		Data(i, knack::AZ) = 0.;
		
		for (size_t j = 0; j < NPARTS; j++) {
			if ( i != j ) {
				partDist = sqrt(	( Data(i, knack::X) - Data(j, knack::X) )*
									( Data(i, knack::X) - Data(j, knack::X) ) +
									( Data(i, knack::Y) - Data(j, knack::Y) )*
									( Data(i, knack::Y) - Data(j, knack::Y) ) +
									( Data(i, knack::Z) - Data(j, knack::Z) )*
									( Data(i, knack::Z) - Data(j, knack::Z) ) );
				partDistPow3 = partDist*partDist*partDist;
				Data(i, knack::AX) -= Data(j, knack::M) * ( Data(i, knack::X) - Data(j, knack::X) )
								/ partDistPow3;
				Data(i, knack::AY) -= Data(j, knack::M) * ( Data(i, knack::Y) - Data(j, knack::Y) )
								/ partDistPow3;
				Data(i, knack::AZ) -= Data(j, knack::M) * ( Data(i, knack::Z) - Data(j, knack::Z) )
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
		
		bfAccLength = sqrt(	BFData(i, knack::AX) * BFData(i, knack::AX)  +
							BFData(i, knack::AY) * BFData(i, knack::AY)  +
							BFData(i, knack::AZ) * BFData(i, knack::AZ) );

		treeAccLength = sqrt(	TreeData(i, knack::AX) * TreeData(i, knack::AX) +
								TreeData(i, knack::AY) * TreeData(i, knack::AY) +
								TreeData(i, knack::AZ) * TreeData(i, knack::AZ) );
							
		accLengthRelErr[i] = ( treeAccLength - bfAccLength ) / bfAccLength;
		std::cout << fabs( accLengthRelErr[i] ) << "\n";
	}
	std::cout << "\n";
	
    MPI::Finalize();
	return EXIT_SUCCESS;
}

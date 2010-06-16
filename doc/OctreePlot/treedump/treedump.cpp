#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/progress.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//#define SINGLEPREC
#define OOSPH_MPI

#ifdef OOSPH_MPI
#include <mpi.h>
#endif

#include "typedefs.h"
#include "octree.h"
#include "nodeproxy.h"
/*typedef NodeProxy	NodeProxyType;
typedef NodeProxy*	NodeProxyPtrType;
typedef NodeProxy&	NodeProxyRefType;*/

//#define NPARTS	0
//#define NPARTS	20
#define NPARTS	50
//#define NPARTS	1000
//#define NPARTS	10000
//#define NPARTS	100000
//#define NPARTS	300000

int main(int argc, char* argv[]) {
	
	MPI::Init(argc, argv);
	
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    const size_t RANK = MPI::COMM_WORLD.Get_rank();
    
    using namespace boost::posix_time;
	ptime TimeStart, TimeStop;
	
	matrixType Data(NPARTS, knack::PSIZE);

	valueType CenterX = 0;
	valueType CenterY = 0;
	for (size_t i = 0; i < NPARTS; i++) {
		Data(i, knack::PID) = i;
		Data(i, knack::X) = ( (valueType)rand() ) / RAND_MAX;
		Data(i, knack::Y) = ( (valueType)rand() ) / RAND_MAX;
		//Data(i, knack::Z) = ( (valueType)rand() ) / RAND_MAX;
		Data(i, knack::Z) = 0.;
		Data(i, knack::M) = 1.;
		CenterX += Data(i, knack::X);
		CenterY += Data(i, knack::Y);
	}
		
	std::vector<knack::NodeProxy> DataProxies;
	DataProxies.resize(NPARTS);
	
	for (size_t i = 0; i < NPARTS; i++) {
		( DataProxies[i] ).setup(&Data, i);
	}

	// start tree context
	TimeStart = microsec_clock::local_time();	
	knack::OctTree BarnesHutTree;

	TimeStop  = microsec_clock::local_time();
	std::cout << "Tree prepare time       " << ( TimeStop - TimeStart ) << "\n";
	
	//boost::progress_display show_progress( NPARTS , std::cout);
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < NPARTS; i++) {
        if ( RANK == 0 ) {
            if ( Data(i, knack::X) < 0.50 ) {
                BarnesHutTree.insertParticle( *(DataProxies[i]),true );
            }
            if ( Data(i, knack::X) > 0.50 && Data(i, knack::X) < 0.75 ) {
                BarnesHutTree.insertParticle( *(DataProxies[i]),false );
            }
        } else {
            if ( Data(i, knack::X) > 0.50 ) {
                BarnesHutTree.insertParticle( *(DataProxies[i]),true );
            }
            if ( Data(i, knack::X) < 0.50 && Data(i, knack::X) > 0.25 ) {
                BarnesHutTree.insertParticle( *(DataProxies[i]),false );
            }
        }
		//BarnesHutTree.insertParticle( *(DataProxies[i]), locality);
        //BarnesHutTree.insertParticle( *(DataProxies[i]), true);
		//++show_progress;
	}
	TimeStop  = microsec_clock::local_time();
	std::cout << "Tree populate time      " << ( TimeStop - TimeStart ) << "\n";
	
    // dump tree
    std::string treeDumpFilename;
        
	TimeStart = microsec_clock::local_time();	
	BarnesHutTree.calcMultipoles();
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Calc. multipoles time   " << ( TimeStop - TimeStart ) << "\n";

    treeDumpFilename = "treedump_";
    treeDumpFilename += boost::lexical_cast<std::string>(RANK+1);
    treeDumpFilename += "_of_";
    treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
    treeDumpFilename += ".dot";
    BarnesHutTree.treeDOTDump(treeDumpFilename);
    
	//show_progress.restart(NPARTS);
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < NPARTS; i++) {
        if ( ( ( RANK == 0 ) && Data(i, knack::X) < 0.50 ) ||
             ( ( RANK == 1 ) && Data(i, knack::X) > 0.50 ) ) {
            BarnesHutTree.calcAcc(*(DataProxies[i]) );
        }
		//++show_progress;
	}
    
    // dump tree
    treeDumpFilename = "treedump_";
    treeDumpFilename += boost::lexical_cast<std::string>(RANK+1);
    treeDumpFilename += "_of_";
    treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
    treeDumpFilename += ".txt";
    BarnesHutTree.treeDump(treeDumpFilename);
    
	TimeStop  = microsec_clock::local_time();
	std::cout << "Gravity calc time       " << ( TimeStop - TimeStart ) << "\n";
	
	MPI::Finalize();
	return EXIT_SUCCESS;
}

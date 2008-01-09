#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdlib>
#include <iostream>
#include <string>

#include <mpi.h>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include <boost/assign/std/vector.hpp>

#include <boost/mpl/vector_c.hpp>

#include "particle.h"


namespace po  = boost::program_options;
namespace mpl = boost::mpl;

#include "simulation_trait.h"
typedef oosph::SimulationTrait<> SimTrait;
typedef SimTrait::value_type value_type;

#include "iomanager.h"
typedef oosph::IOManager<SimTrait> io_type;

#include "memorymanager.h"
typedef oosph::MemoryManager<SimTrait> mem_type;

#include "communicationmanager.h"
typedef oosph::CommunicationManager<SimTrait> com_type;

#include "costzone.h"
typedef mpl::vector_c<size_t, oosph::X> CostZoneIndex;
typedef oosph::CostZone<CostZoneIndex, SimTrait> CostZoneType;

using namespace oosph;
using namespace boost::assign;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <vector>

// tree stuff
#define OOSPH_MPI
#include "octree.h"

int main(int argc, char* argv[])
{
    MPI::Init(argc, argv);

    po::options_description Options("Global Options");
    Options.add_options()("help,h", "Produces this Help blabla...")
    ("input-file,i", po::value<std::string>(), "InputFile");

    po::positional_options_description POD;
    POD.add("input-file",1);

    po::variables_map VMap;
    po::store(po::command_line_parser(argc, argv).options(Options).positional(POD).run(), VMap);
    po::notify(VMap);

    if (VMap.count("help") )
    {
        std::cout << Options << std::endl;
        MPI::Finalize();
        return EXIT_FAILURE;
    }

    if (!VMap.count("input-file") )
    {
        std::cout << Options << std::endl;
        MPI::Finalize();
        return EXIT_FAILURE;
    }

    io_type& IOManager(io_type::Instance() );
    mem_type& MemManager(mem_type::Instance() );
    com_type& ComManager(com_type::Instance() );
    CostZoneType& CostZone(CostZoneType::Instance() );

    SimTrait::matrix_reference Data(MemManager.Data);
    SimTrait::matrix_reference GData(MemManager.GData);
    
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    const size_t RANK = MPI::COMM_WORLD.Get_rank();
    
    std::string InputFileName = VMap["input-file"].as<std::string>();

    Data.resize(Data.size1(), oosph::SIZE);
    GData.resize(GData.size1(), oosph::SIZE);

    IOManager.LoadCDAT(InputFileName);

    std::cout << "Rank " << RANK << ": Number of loaded particles is " << MemManager.Data.size1() << "\n";
    
    CostZone.StartTimer();
    // We're doing some calculations in here...
    sleep(1);
    CostZone.StopTimer();

    ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
    std::cout << "Rank " << RANK <<  ": Local part is " << CostZone.LocalPart()*100 << "% \n";
    
    ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
    std::cout << "Rank " << RANK << ": Composed Ghost Data (" << GData.size1() << ")\n";

    const size_t noParts = Data.size1();
    const size_t noGhosts = GData.size1();
        
    // particles are all distributed now
    using namespace boost::posix_time;
	ptime TimeStart, TimeStop;
    
    std::vector<sphlatch::NodeProxy> partProxies, ghostProxies;
    
    partProxies.resize(noParts);
    for (size_t i = 0; i < noParts; i++) {
        ( partProxies[i] ).setup(&Data, i);
    }
    
    ghostProxies.resize(noGhosts);
    for (size_t i = 0; i < noGhosts; i++) {
        ( ghostProxies[i] ).setup(&GData, i);
    }
    
    TimeStart = microsec_clock::local_time();
    sphlatch::OctTree BarnesHutTree;
    TimeStop  = microsec_clock::local_time();
    std::cerr << "Tree prepare time       " << ( TimeStop - TimeStart ) << "\n";
    
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < noParts; i++) {
        BarnesHutTree.insertParticle( *(partProxies[i]), true);
	}
	for (size_t i = 0; i < noGhosts; i++) {
        //GData(i, M) = 0.; // set to zero check what happens without
                            // ghost gravity
        BarnesHutTree.insertParticle( *(ghostProxies[i]), false);
	}
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Tree populate time      " << ( TimeStop - TimeStart ) << "\n";
    
    std::string treeDumpFilename;

	TimeStart = microsec_clock::local_time();	
	BarnesHutTree.calcMultipoles();    
	TimeStop  = microsec_clock::local_time();
	std::cerr << "Calc. multipoles time   " << ( TimeStop - TimeStart ) << "\n";

    // dump toptree
    treeDumpFilename = "toptree_postsync_";
    treeDumpFilename += boost::lexical_cast<std::string>(RANK+1);
    treeDumpFilename += "_of_";
    treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
    treeDumpFilename += ".txt";
    BarnesHutTree.toptreeDump(treeDumpFilename);
    
    boost::progress_display show_progress( noParts , std::cout);
	TimeStart = microsec_clock::local_time();	
	for (size_t i = 0; i < noParts; i++) {
        BarnesHutTree.calcGravity(*(partProxies[i]) );
		++show_progress;
	}
	TimeStop  = microsec_clock::local_time();
	std::cout << "Gravity calc time       " << ( TimeStop - TimeStart ) << "\n";

    // dump tree;
    treeDumpFilename = "treedump_";
    treeDumpFilename += boost::lexical_cast<std::string>(RANK+1);
    treeDumpFilename += "_of_";
    treeDumpFilename += boost::lexical_cast<std::string>(SIZE);
    treeDumpFilename += ".dot";
    BarnesHutTree.treeDOTDump(treeDumpFilename);

    // save particles
    std::vector<int> outputAttrSet;
    std::string outFilename;
    outFilename += "out";
    //outFilename += "out_noghosts_";
    outFilename += boost::lexical_cast<std::string>(SIZE);
    outFilename += ".cdat";
    outputAttrSet += ID, X, Y, Z, AX, AY, AZ, M;
    IOManager.SaveCDAT(outFilename, outputAttrSet);
        
    MPI::Finalize();
    return EXIT_SUCCESS;
}

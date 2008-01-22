// some defs

#define OOSPH_SINGLE_PRECISION
#define SPHLATCH_SINGLEPREC

//#define OOSPH_LOADBALANCE

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

#include "predictorcorrector.h"
typedef mpl::vector_c<size_t, oosph::X, oosph::VX, oosph::AX, oosph::OX, oosph::OVX, oosph::OAX, oosph::PX, oosph::PVX, oosph::PAX> AccIntIndices;
typedef oosph::SOVecPredictorCorrector<AccIntIndices, SimTrait> AccIntType;

#include "erazer.h"
typedef mpl::vector_c<size_t, oosph::AX, oosph::AY, oosph::AZ> ErazerIndices;
typedef oosph::Erazer<ErazerIndices, SimTrait> ErazerType;

using namespace oosph;
using namespace boost::assign;

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/progress.hpp>
#include <vector>

// tree stuff

#include "bhtree.h"

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
    
    ErazerType& Erazer(ErazerType::Instance() );
    AccIntType AccInt;

    SimTrait::matrix_reference Data(MemManager.Data);
    
    const size_t RANK = MPI::COMM_WORLD.Get_rank();
    
    std::string InputFileName = VMap["input-file"].as<std::string>();

    Data.resize(Data.size1(), oosph::SIZE);

    IOManager.LoadCDAT(InputFileName);
    value_type absTime = 0.;
    
    std::cout << "Rank " << RANK << ": Number of loaded particles is " << MemManager.Data.size1() << "\n";
  
    const size_t noParts = Data.size1();

    // particles are all distributed now
    using namespace boost::posix_time;
	ptime TimeStart, TimeStop;
    
    std::vector<sphlatch::NodeProxy> partProxies;
    
    AccInt.BootStrap();
    
    partProxies.resize(noParts);
    for (size_t i = 0; i < noParts; i++) {
        ( partProxies[i] ).setup(&Data, i);
    }
    
    size_t maxSteps = 150;
    value_type dt = 0.00005;
    
    for (size_t step = 0; step < maxSteps; step++) {
        // predicition context
        { 
            TimeStart = microsec_clock::local_time();
            sphlatch::OctTree BarnesHutTree;
            TimeStop  = microsec_clock::local_time();
            std::cerr << "Tree prepare time       " << ( TimeStop - TimeStart ) << "\n";
    
            TimeStart = microsec_clock::local_time();	
            for (size_t i = 0; i < noParts; i++) {
                BarnesHutTree.insertParticle( *(partProxies[i]), true);
            }
            TimeStop  = microsec_clock::local_time();
            std::cerr << "Tree populate time      " << ( TimeStop - TimeStart ) << "\n";

            TimeStart = microsec_clock::local_time();	
            BarnesHutTree.calcMultipoles();
            TimeStop  = microsec_clock::local_time();
            std::cerr << "Calc. multipoles time   " << ( TimeStop - TimeStart ) << "\n";
    
            boost::progress_display show_progress( noParts , std::cout);
            TimeStart = microsec_clock::local_time();	
            for (size_t i = 0; i < noParts; i++) {
                BarnesHutTree.calcGravity(*(partProxies[i]) );
                ++show_progress;
            }
            TimeStop  = microsec_clock::local_time();
            std::cout << "Gravity calc time       " << ( TimeStop - TimeStart ) << "\n";
        }
        
        // prediction
        AccInt.Predictor(dt);
        absTime += dt;
        MemManager.SaveParameter("TIME", absTime, true);
        std::cerr << "position predicted (" << step << ")\n";
        
        // correction context
        { 
            TimeStart = microsec_clock::local_time();
            sphlatch::OctTree BarnesHutTree;
            for (size_t i = 0; i < noParts; i++) {
                BarnesHutTree.insertParticle( *(partProxies[i]), true);
            }
            BarnesHutTree.calcMultipoles();
            TimeStop  = microsec_clock::local_time();
            std::cerr << "B&H tree build time     " << ( TimeStop - TimeStart ) << "\n";
            
            TimeStart = microsec_clock::local_time();	
            for (size_t i = 0; i < noParts; i++) {
                BarnesHutTree.calcGravity(*(partProxies[i]) );
            }
            TimeStop  = microsec_clock::local_time();
            std::cout << "Gravity calc time       " << ( TimeStop - TimeStart ) << "\n";
        }
        
        // correction
        AccInt.Corrector(dt);
        std::cerr << "position corrected (" << step << ")\n";
        
        // save particles
        if ( ( step % 1 ) == 0 ) {
            std::vector<int> outputAttrSet;
            
            std::string outFilename = "out000000.cdat";
            std::string stepString = boost::lexical_cast<std::string>(step);
            outFilename.replace(outFilename.size() - 5 - stepString.size(),
                stepString.size(), stepString );
            
            outputAttrSet += ID, X, Y, Z, VX, VY, VZ, AX, AY, AZ, M;
            IOManager.SaveCDAT(outFilename, outputAttrSet);
            std::cerr << "saved file " << outFilename << "\n";
        }
        
        Erazer();
    }
    
    MPI::Finalize();
    return EXIT_SUCCESS;
}

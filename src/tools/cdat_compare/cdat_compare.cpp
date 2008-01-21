#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SPH_STANDALONE
#define SPH_STANDALONE
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "simulation_trait.h"
#include "iomanager.h"
#include "particle.h"

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

namespace po = boost::program_options;

typedef oosph::SimulationTrait<> SimTrait;
typedef oosph::IOManager<SimTrait> IOManagerType;
typedef oosph::MemoryManager<SimTrait> MemoryManagerType;
typedef oosph::ParticleVarMap<SimTrait> VarMapType;

typedef SimTrait::matrix_reference matrix_reference;
typedef SimTrait::matrix_type matrix_type;
typedef SimTrait::value_type value_type;

using namespace oosph;
using namespace boost::assign;

int main( int argc, char* argv[] )
{
    IOManagerType& IOManager(IOManagerType::Instance() );
    MemoryManagerType& MemoryManager(MemoryManagerType::Instance() );
    VarMapType& VM(VarMapType::Instance() );
    
    po::options_description Options( "Global Options" );
    Options.add_options()
    ( "ref-file,r", po::value<std::string>(), "reference file" )
    ( "comp-file,c", po::value<std::string>(), "comparison file" );
    
    po::options_description Config( "Config File Options" );
    
    po::positional_options_description p;    
    p.add("ref-file", 1);
    p.add("comp-file", 1);

    po::variables_map VMap;

    po::store( po::command_line_parser( argc, argv ).options( Options ).positional(p).run(), VMap );

    po::notify( VMap );

    if ( VMap.count( "help" ) )
    {
        std::cout << Options << std::endl;
        exit( -1 );
    }

    if ( !VMap.count( "ref-file" ) )
    {
        std::cout << std::endl << "reference file not Specified" << std::endl << std::endl;
        std::cout << Options << std::endl;
        exit( -1 );
    }


    if ( !VMap.count( "comp-file" ) )
    {
        std::cout << std::endl << "comparison file has to be specified" << std::endl << std::endl;
        std::cout << Options << std::endl;
        exit( -1 );
    }

    std::string refFileName( VMap[ "ref-file" ].as<std::string>() );
    std::string compFileName( VMap[ "comp-file" ].as<std::string>() );
	
    std::vector<int> refAttr = IOManager.LoadCDAT(refFileName);    
    matrix_type ref(MemoryManager.Data);
    
    std::vector<int> compAttr = IOManager.LoadCDAT(compFileName);
    matrix_type comp(MemoryManager.Data);
    
    MemoryManager.Data.resize(0,0);
    
    // determine attributes to be compared
    
    std::vector<int> attr;
    for (size_t i = 0; i < compAttr.size(); i++) {
        for (size_t j = 0; j < refAttr.size(); j++) {
            if ( compAttr[i] == refAttr[j] &&
                 compAttr[i] != ID &&
                 compAttr[i] != -1 ) {
                attr += compAttr[i];
                break;
            }
        }
    }
        
    // compare
    const size_t noParts = comp.size1();
    const size_t noRefs  = ref.size1();
    const size_t noAttrs = attr.size();
    
    std::cout << "# " << noParts << " particles in comparison file    "
              << " d_x = x_comp - x_ref \n";
    
    std::cout << "#    id            ";
    for (size_t i = 0; i < noAttrs; i++) {
        std::cout << VM.getName(attr[i]) << "             d_"
                  << VM.getName(attr[i]) << "               ";
    }
    std::cout << "\n";
    
    long int compID, refID;
    
    for (size_t i = 0; i < noParts; i++) {
        compID = static_cast<long int>( comp(i, ID) );
        
        for (size_t j = 0; j < noRefs; j++) {
            refID = static_cast<long int>(  ref(j, ID) );
            if ( compID == refID ) {
                 std::cout << std::setw(7);
                 std::cout << compID << "   ";
                 for (size_t k = 0; k < noAttrs; k++) {
                    std::cout << std::scientific
                              << std::setw(15)
                              << std::setprecision(6)
                              << ref(j, attr[k]) << "   "
                              << comp(i, attr[k]) - ref(j, attr[k])
                              << "   ";
                 }
                 std::cout << "\n";
                 break;
            }
        }
    }
                 
                 
        
   
    // **
    // ** Removing entries outside boundaries
    // **

    //num::matrix<double> Sliced(NumTotParticles, NumAttributes);
    
    
    
    exit( EXIT_SUCCESS );
}


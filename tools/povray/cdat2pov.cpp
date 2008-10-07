/**************************************************************************
*   Copyright (C) 2005 by Pascal Bauer & Andreas Reufer                   *
*   andreas.reufer@space.unibe.ch                                         *
*   pascal.bauer@space.unibe.ch                                           *
*                                                                         *
***************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SPH_STANDALONE
#define SPH_STANDALONE
#endif

#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include "simulation_trait.h"
#include "iomanager.h"
#include "particle.h"

namespace num = boost::numeric::ublas;
namespace po = boost::program_options;

typedef num::matrix<double>::size_type size_type;

typedef oosph::SimulationTrait<> SimTrait;
typedef oosph::IOManager<SimTrait> IOManagerType;
typedef oosph::MemoryManager<SimTrait> MemoryManagerType;
typedef oosph::ParticleVarMap<SimTrait> VarMapType;

typedef SimTrait::value_type value_type;
typedef SimTrait::matrix_reference matrix_reference;
typedef SimTrait::matrix_row matrix_row;

int main( int argc, char* argv[] )
{
    IOManagerType& IOManager(IOManagerType::Instance() );
    MemoryManagerType& MemoryManager(MemoryManagerType::Instance() );
    VarMapType& VM(VarMapType::Instance() );
    
    po::options_description Options( "Global Options" );
    Options.add_options() ( "help", "Produces this Help" )
    ( "input-file,i", po::value<std::string>(), "Input CDAT File" )
    ( "output-file,o", po::value<std::string>(), "Output TXT File" );
            
    po::positional_options_description p;    
    p.add("input-file", 1);
    p.add("output-file", 1);

    po::variables_map VMap;

    po::store( po::command_line_parser( argc, argv ).options( Options ).positional(p).run(), VMap );

    po::notify( VMap );

    if ( VMap.count( "help" ) )
    {
        std::cout << Options << std::endl;
        exit( -1 );
    }

    if ( !VMap.count( "input-file" ) )
    {
        std::cout << std::endl << "Input File not Specified" << std::endl << std::endl;
        std::cout << Options << std::endl;
        exit( -1 );
    }


    if ( !VMap.count( "output-file" ) )
    {
        std::cout << std::endl << "Output File has to be specified" << std::endl << std::endl;
        std::cout << Options << std::endl;
        exit( -1 );
    }
        
    // **
    // ** Start Reading in Data
    // **

    std::string InputFileName( VMap[ "input-file" ].as<std::string>() );
    std::vector<int> ATTRIBUTES = IOManager.LoadCDAT(InputFileName);

    matrix_reference Data(MemoryManager.Data);

    // **
    // ** Start Writing
    // **
    
    std::fstream fout;

    fout.open( VMap[ "output-file" ].as<std::string>().c_str(), std::ios::out );

    if ( !fout ) {
    	std::cerr << "Error Writing " << VMap[ "output-file" ].as<std::string>().c_str() << "\n";
    }
        
    fout << "#version 3.6;\n";
    fout << "global_settings {  assumed_gamma 1.0 }\n";
    fout << "#default{ finish{ ambient 0.1 diffuse 0.9}}\n";
    fout << "#include \"colors.inc\"\n";
    fout << "camera {location <1.2 , 1.2 ,-0.7>\n";
    fout << "        look_at  <0.5 , 0.5 , 0.5>}\n";
    fout << "light_source{<1500,2500,-2500> color White}\n";
    
    for (size_t i = 0; i < Data.size1(); i++) {
        fout << "sphere{<"
                << Data(i, oosph::X)
                << "," << Data(i, oosph::Y)
                << "," << Data(i, oosph::Z)
                << ">, 0.005\n";
        fout << "   texture{\n";
        fout << "       pigment {color rgb<1,0.3,0>} \n";
        fout << "       } \n";
        fout << "} \n";
    }
    
    fout.close();

    std::cout << VMap[ "input-file" ].as<std::string>() << " -> " << VMap[ "output-file" ].as<std::string>() << "\n";
    
    exit( EXIT_SUCCESS );
}


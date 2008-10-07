/**************************************************************************
*   Copyright (C) 2005 by Pascal Bauer & Andreas Reufer                   *
*   andreas.reufer@space.unibe.ch                                         *
*   pascal.bauer@space.unibe.ch                                           *
*                                                                         *
***************************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef OOSPH_STANDALONE
#define OOSPH_STANDALONE
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
    VarMapType& VM(VarMapType::Instance() );
    
    po::options_description Options( "Global Options" );
    Options.add_options() ( "help", "Produces this Help" )
    ( "input-file,i", po::value<std::string>(), "Input CDAT File" )
    ( "output-file,o", po::value<std::string>(), "Output CDAT File" );
    
    po::options_description Config( "Config File Options" );
    
    po::positional_options_description p;    
    p.add("input-file", 1);
    p.add("output-file", 2);

    po::variables_map VMap;

    po::store( po::command_line_parser( argc, argv ).options( Options ).positional(p).run(), VMap );

    po::notify( VMap );

    if ( VMap.count( "help" ) )
    {
        std::cout << Options << std::endl;
        exit( -1 );
    }

    std::string InputFileName;
    if ( !VMap.count( "input-file" ) ) {
        std::cout << "\nInput File not Specified!\n\n" << Options << "\n";
        exit( -1 );
    } else {
	    InputFileName = VMap[ "input-file" ].as<std::string>();
    }

    std::string OutputFileName;
    if ( !VMap.count( "output-file" ) ) {
	    OutputFileName = InputFileName;
    } else {
	    OutputFileName = VMap[ "output-file" ].as<std::string>();
    }

    // Load in any precision
    std::vector<int> ATTRIBUTES = IOManager.LoadCDAT(InputFileName);
    std::cout << InputFileName << " -> " << OutputFileName << "\n";

    for (size_t i = 0; i < ATTRIBUTES.size(); i++) {
	     std::cout << VM.getName(ATTRIBUTES[i]) << " ";
    }
    std::cout << "\n";

    // Save in single precision
    IOManager.SetSinglePrecOut();
    IOManager.SaveCDAT(OutputFileName, ATTRIBUTES);
    
    exit( EXIT_SUCCESS );
}


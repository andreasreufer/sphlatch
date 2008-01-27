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

int main( int argc, char* argv[] )
{
    IOManagerType& IOManager(IOManagerType::Instance() );
    MemoryManagerType& MemoryManager(MemoryManagerType::Instance() );
    VarMapType& VM(VarMapType::Instance() );
 
    using namespace oosph;
    using boost::lexical_cast;

    po::options_description Options( "Global Options" );
    Options.add_options() ( "help", "Produces this Help" )
    ( "input-file,i", po::value<std::string>(), "Input CDAT File" )
    ( "output-file,o", po::value<std::string>(), "Output VTK File" );
    
    po::options_description Config( "Config File Options" );
    
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

    // **
    // ** Checking for coordinates
    // **
    
    using namespace oosph;
    
    std::vector<int> minimums(3);
    minimums[0] = X;
    minimums[1] = Y;
    minimums[2] = Z;

    for (size_t i = 0; i < minimums.size(); ++i) {
    	int sizepreerase =  ATTRIBUTES.size();
    	ATTRIBUTES.erase(remove(ATTRIBUTES.begin(), ATTRIBUTES.end(), minimums[i]),ATTRIBUTES.end());
    	int sizeposterase =  ATTRIBUTES.size();
	if (sizepreerase - sizeposterase == 0) {
		std::cerr << "Variable \"" << VM.getName(minimums[i]) << "\" not found in inputfile, aborting!\n";
		exit( EXIT_FAILURE );
		};
	}
  
    // **
    // ** Start Writing
    // **
    
    std::fstream fout;
    std::string writtenattribs;
    int numtotparticles = MemoryManager.Data.size1();
    int first, second, third;

    fout.open( VMap[ "output-file" ].as<std::string>().c_str(), std::ios::out | std::ios::ate);

    if ( !fout ) {
    	std::cerr << "Error Writing " << VMap[ "output-file" ].as<std::string>().c_str() << "\n";
    }

    fout << "# vtk DataFile Version 2.0\n";
    fout << MemoryManager.Name << "\n";
    fout << "ASCII\n";
    fout << "DATASET UNSTRUCTURED_GRID\n";
    fout << "POINTS " << lexical_cast<std::string>(numtotparticles) << " double\n";
    fout << std::scientific;

    first = minimums[0];
    second = minimums[1];
    third = minimums[2];

    for (int i = 0; i < numtotparticles; ++i) {
	fout << MemoryManager.Data(i, first) << "\t" << MemoryManager.Data(i, second) << "\t" << MemoryManager.Data(i, third) << "\n";
	}
  
    fout << "POINT_DATA " << lexical_cast<std::string>(numtotparticles) << "\n";
    std::vector<int> next(27);
    next[ID] = -1;
    next[X] = Y;
    next[Y] = Z;
    next[Z] = X;
    next[VX] = VY;
    next[VY] = VZ;
    next[VZ] = VX;
    next[AX] = AY;
    next[AY] = AZ;
    next[AZ] = AX;
    next[M] = -1;
    next[H] = -1;
    next[DHDT] = -1;
    next[RHO] = -1;
    next[E] = -1;
    next[P] = -1;
    next[POW] = -1;
    next[DIV_V] = -1;
    next[ROTX_V] = ROTY_V;
    next[ROTY_V] = ROTZ_V;
    next[ROTZ_V] = ROTX_V;
    next[Q] = -1;
    next[GRAVEPS] = -1;
    next[MAXAVMU] = -1;
    next[ALPHA] = -1;
    next[DALPHADT] = -1;
    next[NONEIGH] = -1;

    while(!ATTRIBUTES.empty()) {
    	first = ATTRIBUTES.front();
    	ATTRIBUTES.erase(remove(ATTRIBUTES.begin(), ATTRIBUTES.end(), first),ATTRIBUTES.end());
    	second = third = -1;
    	for (size_t i = 0; i < ATTRIBUTES.size(); ++i) {
    		if (ATTRIBUTES[i] == next[first]) {
    			second = next[first];
    			ATTRIBUTES.erase(remove(ATTRIBUTES.begin(), ATTRIBUTES.end(), second),ATTRIBUTES.end());
    			for (size_t j = 0; j < ATTRIBUTES.size(); ++j) {
    				if (ATTRIBUTES[j] == next[second]) {
    				third = next[second];
    				ATTRIBUTES.erase(remove(ATTRIBUTES.begin(), ATTRIBUTES.end(), third),ATTRIBUTES.end());
    				}
    			}
    		}
    	}
	if (second != -1) {
		// Vector output
		writtenattribs += VM.getName(first);
		writtenattribs += VM.getName(second);
		writtenattribs += VM.getName(third);
		writtenattribs += " ";
		fout << "VECTORS " << VM.getName(first) << VM.getName(second) << VM.getName(third) << " double\n";
    			for (int i = 0; i < numtotparticles; ++i) {
				fout << MemoryManager.Data(i, first) << "\t" << MemoryManager.Data(i, second) << "\t" << MemoryManager.Data(i, third) << "\n";
			}

		} else {
		// Scalar output
		writtenattribs += VM.getName(first);
		writtenattribs += " ";
		fout << "SCALARS " << VM.getName(first) << " double\n";
		fout << "LOOKUP_TABLE default\n";
    			for (int i = 0; i < numtotparticles; ++i) {
				fout << MemoryManager.Data(i, first) << "\n";
				}
		}
    }
			

    fout.close();
   
    std::cout << VMap[ "input-file" ].as<std::string>() << " -> " << VMap[ "output-file" ].as<std::string>() << "\n";
    std::cout << "Written attributes: " << writtenattribs << "\n";

    exit( EXIT_SUCCESS );
}


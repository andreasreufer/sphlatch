/**************************************************************************
*   2006 by Andreas Reufer                                                *
*   andreas.reufer@space.unibe.ch                                         *
*                                                                         *
***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef OOSPH_STANDALONE
#define OOSPH_STANDALONE
#endif

#include <string>
#include <vector>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/progress.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <hdf5.h>

#include "simulation_trait.h"
#include "iomanager.h"
#include "particle.h"

namespace po = boost::program_options;
namespace mpl = boost::mpl;

typedef oosph::SimulationTrait<> SimTrait;
typedef oosph::IOManager<SimTrait> IOManagerType;
typedef oosph::MemoryManager<SimTrait> MemoryManagerType;
typedef oosph::ParticleVarMap<SimTrait> VarMapType;

typedef SimTrait::value_type value_type;
typedef SimTrait::vector_type vector_type;
typedef SimTrait::matrix_type matrix_type;
typedef SimTrait::matrix_reference matrix_reference;

typedef SimTrait::matrix_column matrix_column;
typedef SimTrait::matrix_row matrix_row;

typedef SimTrait::matrix_row_range matrix_row_range;

using namespace oosph;

int main( int argc, char* argv[] )
{
    IOManagerType& IOManager(IOManagerType::Instance() );
    MemoryManagerType& MemoryManager(MemoryManagerType::Instance() );
    VarMapType& VM(VarMapType::Instance()); 
    
    std::string InputFileName, OutputFileName;

    po::options_description Options( "Global Options" );
    Options.add_options() ( "help", "Produces this Help" )
    ( "input-file,i", po::value<std::string>(),  "Input CDAT File" )
    ( "output-file,o", po::value<std::string>(), "Output HDF5 Data File" )
    ;
    
    po::options_description Config( "Config File Options" );
    
    po::positional_options_description p;    
    p.add("input-file", 1);
    p.add("output-file", 1);

    po::variables_map VMap;

    po::store( po::command_line_parser( argc, argv ).options( Options ).run(), VMap );

    po::notify( VMap );

    if ( VMap.count( "help" ) ) {
        std::cout << Options << std::endl;
        exit( -1 );
    }

    std::string InputFilename, OutputFilename;
    
    if ( !VMap.count( "input-file" ) ) {
	    std::cout << "\nInput File not Specified\n\n";
	    std::cout << Options << "\n";
        exit( -1 );
    } else {
	    InputFilename = VMap[ "input-file" ].as<std::string>();
    }

    if ( !VMap.count( "output-file" ) ) {
	    std::cout << "\nOutput File not Specified\n\n";
	    std::cout << Options << "\n";
        exit( -1 );
    } else {
	    OutputFilename = VMap[ "output-file" ].as<std::string>();
    }

    matrix_reference Data(MemoryManager.Data);


    // Prepare HDF5 File

    hid_t file_id = H5Fcreate ( OutputFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[2] = { Data.size1(), 1 };


    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

    for (size_t i = 0; i < ATTRIBUTES.size(); i++) {
	    std::string VarName = VM.getName(ATTRIBUTES[i]);

	    hid_t var_dataset_id = H5Dcreate(file_id, VarName.c_str(), OutputDatatype, dataspace_id, H5P_DEFAULT);
	    H5Dclose( var_dataset_id );
    }

    std::cout << "\n";
    
    H5Fclose( file_id );

    exit( EXIT_SUCCESS );
}


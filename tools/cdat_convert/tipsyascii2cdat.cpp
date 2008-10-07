/**************************************************************************
*   Copyright (C) 2006 by Andreas Reufer                                  *
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
#include <boost/lexical_cast.hpp>
#include <boost/assign/std/vector.hpp>

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
    
    po::options_description Options( "Global Options" );
    Options.add_options() ( "help", "Produces this Help" )
    ( "input-file,i", po::value<std::string>(),   "Input Tipsy ASCII File" )
    ( "output-file,o", po::value<std::string>(),  "Output CDAT File" )
    ( "dataset-name,n", po::value<std::string>(), "Name of CDAT dataset         (default:   <none>)" )
    ( "energy-factor,f", po::value<double>(),     "Conversion factor for specific energy u \nu = fact / ( 1 - gamma ) * T (default: 8.58E-3 )" )
    ( "star-energy,U", po::value<double>(),       "Energy for stars             (default:       0 )" )
    ( "gas-gamma,g", po::value<double>(),         "Gamma for EOS                (default:     5/3 )" );
    
    po::options_description Config( "Config File Options" );
    
    po::positional_options_description p;    
    p.add("input-file", 1);
    p.add("output-file", 1);
    p.add("dataset-name", 1);
    p.add("energy-factor", 1);
    p.add("star-energy", 1);
    p.add("gas-gamma", 1);

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

    if ( !VMap.count( "dataset-name") )
    {
	MemoryManager.Name = "<none>";
    } else {
	MemoryManager.Name =  VMap[ "dataset-name" ].as<std::string>();
    }

    double StarEnergy;

    if ( !VMap.count( "star-energy") )
    {
	StarEnergy = 0;
    } else {
	StarEnergy =  VMap[ "star-energy" ].as<double>();
    }
    
    double GasGamma;

    if ( !VMap.count( "gas-gamma") )
    {
	GasGamma = 5./3.;
    } else {
	GasGamma =  VMap[ "gas-gamma" ].as<double>();
    }

    double EnergyConversion;

    if ( !VMap.count( "energy-factor") )
    {
	EnergyConversion = 8.5861E-3 / ( GasGamma - 1 );
    } else {
	EnergyConversion =  VMap[ "energy-factor" ].as<double>() / ( GasGamma - 1 );
    }

    // **
    // ** Start Reading in Data
    // **
    
    std::cout << VMap[ "input-file" ].as<std::string>() << " -> " << VMap[ "output-file" ].as<std::string>() << "\n\n";
    
    std::string InputFileName( VMap[ "input-file" ].as<std::string>() );
    std::fstream fin;

    fin.open( VMap[ "input-file" ].as<std::string>().c_str(), std::ios::in );

    if ( !fin ) {
    	std::cerr << "Error Writing " << VMap[ "input-file" ].as<std::string>().c_str() << "\n";
    }

    std::string ReadLine;
    size_t ntotal, ngas, nstar, ndark;
    
    using boost::lexical_cast;
    using namespace boost::assign;
    using namespace oosph;
    
    fin >> ReadLine; ntotal = lexical_cast<size_t>(ReadLine);
    fin >> ReadLine; ngas = lexical_cast<size_t>(ReadLine);
    fin >> ReadLine; nstar = lexical_cast<size_t>(ReadLine);
    ndark = ntotal - ngas - nstar;

    std::cout << "Number of gas particles:         " << ngas << "\n"
	      << "Number of dark matter particles: " << ndark << "\n"
	      << "Number of stars:                 " << nstar << "\n\n";

    MemoryManager.Data.resize(ntotal, SIZE);

    fin >> ReadLine;
    
    fin >> ReadLine;
    MemoryManager.SaveParameter("TIME", lexical_cast<double>(ReadLine), true);
    MemoryManager.SaveParameter("GAMMA", GasGamma, true);

    std::vector<int> TOTPOSVEL;
    TOTPOSVEL += M, X, Y, Z, VX, VY, VZ;

    for (size_t i = 0; i < TOTPOSVEL.size(); i++) {
	    size_t ATTR = TOTPOSVEL[i];
	    for (size_t j = 0; j < ntotal; j++) {
		  fin >> ReadLine; // Reading mass, position and velocity for ALL particles SaveParameter("GAMMA", gamma, true);
		  MemoryManager.Data(j, ATTR) = lexical_cast<double>(ReadLine);
	    }
    }

    for (size_t j = ngas; j < ntotal; j++) {
	    fin >> ReadLine;	// Reading gravitation softening length for STARS and DARK MATTER
	    MemoryManager.Data(j, GRAVEPS) = lexical_cast<double>(ReadLine);
	    std::cout << "Star or dark matter particle with ID " << j + 1 << "\n";
    }

    /*std::vector<int> GAS;
    GAS += RHO, E, GRAVEPS;

    for (size_t i = 0; i < GAS.size(); i++) {
	    size_t ATTR = GAS[i];
	    for (size_t j = 0; j < ngas; j++) {
		  fin >> ReadLine; // Reading density, temperature and gravitational smoothing length for GAS
		  MemoryManager.Data(j, ATTR) = lexical_cast<double>(ReadLine);
	    }
    }*/
    
    for (size_t j = 0; j < ngas; j++) {
	    fin >> ReadLine; // Reading density
	    MemoryManager.Data(j, RHO) = lexical_cast<double>(ReadLine);
    }
    
    for (size_t j = 0; j < ngas; j++) {
	    fin >> ReadLine; // Reading temperature
	    MemoryManager.Data(j, E) = EnergyConversion * lexical_cast<double>(ReadLine);
    }

    for (size_t j = 0; j < ngas; j++) {
	    fin >> ReadLine; // Reading gravtitation epsilon
	    MemoryManager.Data(j, GRAVEPS) = lexical_cast<double>(ReadLine);
    }

    for (size_t j = 0; j < ntotal; j++) {
	    MemoryManager.Data(j, ID) = j + 1; // Set missing things for ALL particles
    }

    for (size_t j = 0; j < ngas; j++) {
	    // Set missing things for GAS
    }

    for (size_t j = ngas; j < ntotal; j++) {
	    MemoryManager.Data(j, RHO) = 0; // Set missing things for STARS and DARK MATTER
	    MemoryManager.Data(j, E) = StarEnergy;
    }

    fin.close();

    std::vector<int> OutputAttrSet;
    OutputAttrSet += ID, X, Y, Z, VX, VY, VZ, M, E, GRAVEPS;

    std::string OutputFileName( VMap[ "output-file" ].as<std::string>() );

    IOManager.SaveCDAT(OutputFileName, OutputAttrSet);

    exit( EXIT_SUCCESS );
}


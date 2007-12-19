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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/assign/std/vector.hpp>

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

using namespace std;
using namespace oosph;
using namespace boost::assign;

int main(int argc, char* argv[])
{
    // Get the Command line Arguments

    po::options_description Options("Global Options");

    Options.add_options()
    ("help,h", "Produce this help Blabla...")
    ("output-file,o", po::value<std::string>(), "Output File Name")
    ("number-of-particles,n", po::value<size_t>(), "Number of Particles")
    ("neighbours-per-particle,b", po::value<size_t>(), "Neighbours per Particles")
    ("radius,r", po::value<double>(), "Radius of the Initial Sphere");

    po::positional_options_description POD;
    POD.add("output-file", 1);
    POD.add("number-of-particles", 1);

    po::variables_map VMap;
    po::store(po::command_line_parser(argc, argv).options(Options).positional(POD).run(), VMap);
    po::notify(VMap);

    // Check The Arguments

    if (VMap.count("help") )
    {
        std::cout << Options << std::endl;
        return EXIT_FAILURE;
    }

    std::string OutputFileName;

    if (VMap.count("output-file") )
    {
        OutputFileName = VMap["output-file"].as<std::string>();
    }
    else
    {
        std::cout << Options << std::endl;
        return EXIT_FAILURE;
    }

    size_t NOP;

    if (VMap.count("number-of-particles") )
    {
        NOP = VMap["number-of-particles"].as<size_t>();
    }
    else
    {
        std::cout << Options << std::endl;
        std::cout << "Number of Particles must be specified" << std::endl;
        return EXIT_FAILURE;
    }

    size_t NPP;

    if (VMap.count("neighbours-per-particle") )
    {
        NPP = VMap["neighbours-per-particle"].as<size_t>();
    }
    else
    {
        NPP = 50;
    }

    double Radius;

    if (VMap.count("radius") )
    {
        Radius = VMap["radius"].as<double>();
    }
    else
    {
        Radius = 1.;
    }

    IOManagerType& IOManager(IOManagerType::Instance() );
    MemoryManagerType& MemoryManager(MemoryManagerType::Instance() );
    
    MemoryManager.Data.resize(NOP, 20);

    
    const double Volume = 4. / 3. * M_PI * pow(Radius, 3.);
    const double ParticleRadius = Radius / cbrt(1.663 * NOP) * 5;

    const double cubesize = 1;
    const double smoothlength = pow( ((double)NPP/NOP), 1./3. ) * cubesize;

    for (size_t i = 0; i < MemoryManager.Data.size1(); i++) {
    	MemoryManager.Data(i, ID) = i + 1;
    	MemoryManager.Data(i, X) = cubesize*((double)rand() / RAND_MAX);
    	MemoryManager.Data(i, Y) = cubesize*((double)rand() / RAND_MAX);
    	MemoryManager.Data(i, Z) = cubesize*((double)rand() / RAND_MAX);
	MemoryManager.Data(i, VX) = 0.;
	MemoryManager.Data(i, VY) = 0.;
	MemoryManager.Data(i, VZ) = 0.;
    	MemoryManager.Data(i, M) = 1.;
    	MemoryManager.Data(i, E) = 10.*((double)rand() / RAND_MAX) + 1;
    	MemoryManager.Data(i, H) = smoothlength;
    }
    
    MemoryManager.Name = "RandomParticles in a Cube";

    std::vector<int> SaveAttributes;
    SaveAttributes += ID, X, Y, Z, VX, VY, VZ, H, M;
    IOManager.SaveCDAT(OutputFileName, SaveAttributes);

    std::cout << "Wrote " << NOP << " with a smoothing length of " << smoothlength << "\n";

    return EXIT_SUCCESS;
}


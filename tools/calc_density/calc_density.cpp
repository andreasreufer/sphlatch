#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define OOSPH_SINGLE_PRECISION

#include <mpi.h>

#include <cstdlib>
#include <iostream>
#include <string>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>

#include <boost/assign/std/vector.hpp>

#include <boost/mpl/vector_c.hpp>

#include "particle.h"

#include "executor.h"

namespace po  = boost::program_options;
namespace mpl = boost::mpl;

#include "simulation_trait.h"
typedef oosph::SimulationTrait<> SimTrait;

#include "iomanager.h"
typedef oosph::IOManager<SimTrait> io_type;

#include "memorymanager.h"
typedef oosph::MemoryManager<SimTrait> mem_type;

#include "communicationmanager.h"
typedef oosph::CommunicationManager<SimTrait> com_type;

#include "costzone.h"
typedef mpl::vector_c<size_t, oosph::X> CostZoneIndex;
typedef oosph::CostZone<CostZoneIndex, SimTrait> CostZoneType;

#include "rankspace.h"
typedef mpl::vector_c<size_t, oosph::X, oosph::H, oosph::NONEIGH> rankspace_index;
typedef oosph::RankSpace<rankspace_index, SimTrait> neighbour_type;

/** Define the Kernel */
typedef mpl::vector_c<size_t, oosph::X, oosph::H> KernelIndices;
/*#include "cubicspline.h"
typedef oosph::CubicSpline<KernelIndices, SimTrait> kernel_type;*/
#include "symcubicspline.h"
typedef oosph::SymCubicSpline<KernelIndices, SimTrait> kernel_type;


/** Define the Artificial Viscosity */
#include "monaghan.h"
typedef mpl::vector_c<size_t, oosph::X, oosph::VX, oosph::RHO, oosph::H, oosph::P, oosph::MAXAVMU> AVIndices;
typedef oosph::Monaghan<AVIndices, SimTrait> av_type;

/** Define the Equation of State */
#include "polytropicidealgas.h"
typedef mpl::vector_c<size_t, oosph::RHO, oosph::E, oosph::P> EOSIndices;
typedef oosph::PolytropicIdealGas<EOSIndices, SimTrait> eos_type;

#include "sph_trait.h"
typedef oosph::SPH_Trait<kernel_type, av_type, eos_type, neighbour_type> SPHTrait;

#include "densityfunction.h"
typedef mpl::vector_c<size_t, oosph::M, oosph::RHO> DensityIndices;
typedef oosph::DensityFunction<DensityIndices, SPHTrait> density_type;

typedef mpl::vector<density_type*> DensQueue;
typedef oosph::Executor<DensQueue, SPHTrait> DensExecutorType;

typedef SimTrait::value_type value_type;
typedef SimTrait::index_vector_type index_vector_type;

using namespace oosph;
using namespace boost::assign;

int main(int argc, char* argv[])
{
    MPI::Init(argc, argv);

    po::options_description Options("Global Options");
    Options.add_options()("help,h", "Produces this Help")
    ("input-file,i", po::value<std::string>(), "InputFile")
    ("output-file,o", po::value<std::string>(), "OutputFile");

    po::positional_options_description POD;
    POD.add("input-file",1);
    POD.add("output-file",1);

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

    std::string OutputFile;

    if (!VMap.count("output-file") )
    {
        OutputFile = "Data";
    }
    else
    {
        OutputFile = VMap["output-file"].as<std::string>();
    }


    io_type& IOManager(io_type::Instance() );
    mem_type& MemManager(mem_type::Instance() );
    com_type& ComManager(com_type::Instance() );

    SimTrait::matrix_reference Data(MemManager.Data);
    SimTrait::matrix_reference GData(MemManager.GData);

    const size_t RANK = MPI::COMM_WORLD.Get_rank();
    const size_t SIZE = MPI::COMM_WORLD.Get_size();
    
    std::string InputFileName = VMap["input-file"].as<std::string>();

    Data.resize(Data.size1(), oosph::SIZE);
    // Ghosts don't need any integration variables
    GData.resize(GData.size1(), oosph::OX);

    std::vector<int> InputAttrSet = IOManager.LoadCDAT(InputFileName);
    std::cout << "Rank " << RANK << ": Number of loaded particles is " << MemManager.Data.size1() << "\n";

    // Instaniate classes which may depend on input data parameters
    //eos_type& EOS(eos_type::Instance() );
    CostZoneType& CostZone(CostZoneType::Instance() );
    neighbour_type& Neighbours(neighbour_type::Instance() );
    DensExecutorType& DensExecutor(DensExecutorType::Instance() );
    
    ComManager.Exchange(Data, CostZone.CreateDomainIndexVector(), Data);
    std::cout << "Rank " << RANK << ": composed Local Data (" << CostZone.LocalPart()*100 << "%)\n";

    ComManager.Exchange(Data, CostZone.CreateDomainGhostIndexVector(), GData);
    std::cout << "Rank " << RANK << ": Composed Ghost Data (" << MemManager.GData.size1() << ")\n";

    std::cout << "Rank " << RANK << ": neighbour search\n";
    Neighbours();
    std::cout << "Rank " << RANK << ": density, velocity divergence, vorticity\n";
    DensExecutor();
    
	std::string OutputName = OutputFile + "000000.cdat";
	InputAttrSet += RHO;
	IOManager.SaveCDAT(OutputFile, InputAttrSet);
    
    MPI::Finalize();
    return EXIT_SUCCESS;
}

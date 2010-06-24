#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP
//#include <mpi.h>
//#define SPHLATCH_MPI

#define SPHLATCH_HDF5

#include "typedefs.h"
typedef sphlatch::fType   fType;

#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart
{
public:
   // register all variables, which should be loaded from the input file
   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(id, "id"));
      
      vars.push_back(storeVar(h, "h"));

      return(vars);
   }
   
   // register all variables, which should be loaded from the input file
   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      return(vars);
   }
};


class ghost :
   public sphlatch::treeGhost,
   public sphlatch::movingGhost,
   public sphlatch::SPHfluidGhost { };

typedef particle                       partT;
typedef ghost                          ghstT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;


int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   partSetT particles;

   double start;
   start = omp_get_wtime();
   particles.loadHDF5("twophase.h5part");
   std::cout << "particles loaded " << omp_get_wtime() - start << "s\n";
   
   
   start = omp_get_wtime();
   std::cout << particles.step << "\n";
   particles.saveHDF5("out.h5part");
   std::cout << "particles stored " << omp_get_wtime() - start << "s\n";

   particles.step = 42;
   std::cout << particles.step << "\n";
   particles.saveHDF5("multistep.h5part");
   
   particles.step = 43;
   std::cout << particles.step << "\n";
   particles.saveHDF5("multistep.h5part");


#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

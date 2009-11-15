#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP
//#include <mpi.h>
//#define SPHLATCH_MPI

#include "typedefs.h"
typedef sphlatch::fType   fType;

/*#include "bhtree.cpp"
   typedef sphlatch::BHTree     treeT;

   typedef sphlatch::pnodT      pnodT;
   typedef sphlatch::pnodPtrT   pnodPtrT;

   typedef sphlatch::nodeT      nodeT;
   typedef sphlatch::nodePtrT   nodePtrT;

   typedef sphlatch::gcllT      gcllT;
   typedef sphlatch::gcllPtrT   gcllPtrT;

   typedef sphlatch::czllT      czllT;
   typedef sphlatch::czllPtrT   czllPtrT;*/


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

      return(vars);
   }
};


class ghost :
   public sphlatch::treeGhost,
   public sphlatch::movingGhost,
   public sphlatch::SPHfluidGhost { };

typedef particle                       partT;
typedef ghost                          ghstT;

#define SPHLATCH_HDF5
#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>   partSetT;

/*#include "bhtree_treedump.cpp"
   typedef sphlatch::BHTreeDump                   dumpT;

 #include "bhtree_worker_grav.cpp"
   typedef sphlatch::fixThetaMAC                  macT;
   typedef sphlatch::GravityWorker<macT, partT>   gravT;*/

//#include "log_manager.h"
//typedef sphlatch::LogManager                   logT;

//#include "bhtree_worker_sphsum.cpp"


int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif
   //treeT& Tree(treeT::instance());

   const fType  costMin = 1.0e4, costMax = 1.5e4;
   const size_t noParts = 1000000;

   //const size_t noParts = 300000;

   //const fType costMin = 10., costMax = 15.;
   //const size_t noParts = 100;

   partSetT particles;
   particles.resize(noParts);

   particles.loadDump("test.hdf5");

   std::cout << particles.step << "\n";

   for (size_t i = 0; i < noParts; i++)
   {
      particles[i].pos[0] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[1] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[2] = static_cast<fType>(rand()) / RAND_MAX;

      particles[i].m = 1.;

      particles[i].id   = i;
      particles[i].cost = 1.;
   }

   double start;
   start = omp_get_wtime();
   for (size_t i = 0; i < noParts; i++)
   {
      //Tree.insertPart(particles[i]);
   }
   std::cout << "particles insert " << omp_get_wtime() - start << "s\n";


#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

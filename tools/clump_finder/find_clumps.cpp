#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::cType     cType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

const fType finf = sphlatch::fTypeInf;

#include "bhtree.cpp"
typedef sphlatch::BHTree    treeT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
#include "clump_particle.h"
class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::clumpPart,
   public sphlatch::IOPart
{
public:

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(rho, "rho"));
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(clumpid, "clumpid"));
#ifdef CLUMPFIND_FOF
#else
      vars.push_back(storeVar(pot, "pot"));
#endif
      return(vars);
   }
};

typedef particle                               partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "clump_finder.cpp"
typedef sphlatch::Clumps<partT>                clumpsT;

#include "bhtree.cpp"
#include "bhtree_worker_grav.cpp"
typedef sphlatch::BHTree                       treeT;
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

// particles are global
partSetT parts;
clumpsT  clumps;

int main(int argc, char* argv[])
{
   if (not ((argc == 5) || (argc == 6)))
   {
      std::cerr <<
      "usage: clump_finder <inputdump> <clumpsfile> <minrho> <hmult> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string clFilename = argv[2];

   if (argc == 4)
   {
      std::istringstream threadStr(argv[3]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }


   // load the particles
   parts.loadHDF5(inFilename);

   const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(_parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   Tree.update(0.8, 1.2);

#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcAcc(CZbottomLoc[i]);

   clumps.getClumpsPot(parts);
   parts.saveHDF5(inFilename);


   std::cerr << clumps.getNop() << " clump(s) found!\n";
   clumps.doublePrecOut();
   clumps.saveHDF5(clFilename);

   return(0);
}

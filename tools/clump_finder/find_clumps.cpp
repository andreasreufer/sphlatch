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
      vars.push_back(storeVar(pot, "pot"));
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
   if (not ((argc == 4) || (argc == 5)))
   {
      std::cerr <<
      "usage: clump_finder <inputdump> <clumpsfile> <minmass> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string clFilename = argv[2];
   
   std::istringstream mMassStr(argv[3]);
   fType minMass;
   mMassStr >> minMass;

   if (argc == 5)
   {
      std::istringstream threadStr(argv[4]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   // load the particles
   parts.loadHDF5(inFilename);

   const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;
   const fType  G         = parts.attributes["gravconst"];

   treeT& Tree(treeT::instance());
   gravT  gravWorker(&Tree, G);

   Tree.setExtent(parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   Tree.update(0.8, 1.2);

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcPot(CZbottomLoc[i]);

   clumps.getClumps(parts, minMass);
   std::cerr << clumps.getNop()-1 << " clump(s) found!\n";

   parts.saveHDF5(inFilename);
   clumps.doublePrecOut();
   clumps.saveHDF5(clFilename);

   return(0);
}

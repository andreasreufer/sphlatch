#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH
#define SPHLATCH_GRAVITY_POTENTIAL

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::cType     cType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

const fType finf = sphlatch::fTypeInf;

#include "bhtree.cpp"
typedef sphlatch::BHTree   treeT;

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
      vars.push_back(storeVar(pot, "pot"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;
      
      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      
      vars.push_back(storeVar(clumpid, "clumpid"));
      vars.push_back(storeVar(orbit, "orbit"));
      
      vars.push_back(storeVar(rho, "rho"));
      return vars;
   }
};

typedef particle   partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "clump_finder.cpp"
typedef sphlatch::Clumps<partT> clumpsT;

// particles are global
partSetT parts;
clumpsT clumps;

int main(int argc, char* argv[])
{
   if (not ((argc == 2) || (argc == 3)))
   {
      std::cerr <<
      "usage: clump_finder <inputdump> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];

   if (argc == 3)
   {
      std::istringstream threadStr(argv[2]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   // load the particles
   parts.loadHDF5(inFilename);
   const size_t nop = parts.getNop();

   fType pMinMass = 0.;
   for (size_t i = 0; i < nop; i++)
      pMinMass = parts[i].m > pMinMass ? parts[i].m : pMinMass;
      const fType cMinMass = 10. * pMinMass;

   clumps.getClumps(parts, cMinMass);

   clumps[1].getCentralBodyOrbits(0.5, parts.attributes["gravconst"]);
   parts.saveHDF5(inFilename);

   clumps.doublePrecOut();
   clumps.saveHDF5("clump.h5part");
   

   return(0);
}

#include <iostream>
#include <sstream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"

typedef sphlatch::attrMT attrMT;

///
/// define the particle we are using
///
#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"
#include "integrator_predcorr.cpp"

#ifdef SPHLATCH_FIND_CLUMPS
#include "clump_particle.h"
#endif

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::IOPart
#ifdef SPHLATCH_ANEOS
   , public sphlatch::ANEOSPart
#endif
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

      vars.push_back(storeVar(u, "u"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
#ifdef SPHLATCH_GRAVITY_EPSSMOOTHING
      vars.push_back(storeVar(eps, "eps"));
#endif
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

      vars.push_back(storeVar(u, "u"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
#ifdef SPHLATCH_GRAVITY_EPSSMOOTHING
      vars.push_back(storeVar(eps, "eps"));
#endif
      return(vars);
   }
};

typedef particle   partT;

class ghost :
   public sphlatch::treeGhost,
   public sphlatch::movingGhost
{ };

typedef ghost                                  ghstT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;



int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not (argc == 4) )
   {
      std::cerr << "usage: h5part_combine__ <dumpA> <dumpB> <outputDump>\n";
      return(1);
   }

   partSetT partsA, partsB, partsO;

   partsA.loadHDF5(argv[1]);
   partsB.loadHDF5(argv[2]);

   const size_t noA = partsA.getNop();
   const size_t noB = partsB.getNop();


   partsO.resize(noA + noB);

   for (size_t i = 0; i < noA; i++)
     partsO[i] = partsA[i];
   
   attrMT::const_iterator aItr;
   for (aItr  = partsA.attributes.begin();
        aItr != partsA.attributes.end();
        aItr++)
     partsO.attributes.insert(*aItr);

   for (size_t i = 0; i < noB; i++)
   {
     partsO[i+noA] = partsB[i];
   }

   for (aItr  = partsB.attributes.begin();
        aItr != partsB.attributes.end();
        aItr++)
     partsO.attributes.insert(*aItr);
   
   partsO.step = partsA.step;
   partsO.saveHDF5(argv[3]);

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

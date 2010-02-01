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

#define SPHLATCH_LOGGER
#include "logger.cpp"
typedef sphlatch::Logger   logT;

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

      vars.push_back(storeVar(u, "u"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;
      return vars;
   }
};

typedef particle   partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "clumps.cpp"
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
      std::istringstream threadStr(argv[4]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   logT& Logger(logT::instance());

   ///
   /// log program compilation time
   ///
   Logger.stream << "executable compiled from " << __FILE__
                 << " on " << __DATE__
                 << " at " << __TIME__ << "\n\n"
                 << "    features: \n"
                 << "     basic SPH\n";
   Logger.flushStream();

   Logger.stream << "working on " << omp_get_num_threads() << " threads";
   Logger.flushStream();

   // load the particles
   parts.loadHDF5(inFilename);
   Logger << "loaded particles";

   fType& time(parts.attributes["time"]);
   cType& step(parts.step);
   treeT& Tree(treeT::instance());

   const size_t nop = parts.getNop();

   // first bootstrapping step 

   //clumps.clear();
   clumps.getClumps(parts, 1.0);

   /*clumpsT::const_iterator cItr;
   for (cItr = clumps.begin(); cItr != clumps.end(); cItr++)
     std::cout << (*cItr).m << "\n";*/

   std::cout << clumps.getNop() << "\n";

   clumps.doublePrecOut();
   clumps.saveHDF5("clump.h5part");

   return(0);
}

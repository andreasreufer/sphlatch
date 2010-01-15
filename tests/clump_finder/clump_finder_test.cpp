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

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
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

      vars.push_back(storeVar(u, "u"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
#endif
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(acc, "acc"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(noneigh, "noneigh"));

      vars.push_back(storeVar(u, "u"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(T, "T"));
      vars.push_back(storeVar(phase, "phase"));
#endif
      
      vars.push_back(storeVar(cost, "cost"));
      return(vars);
   }
};

typedef particle   partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "sph_algorithms.cpp"
#include "sph_kernels.cpp"
#include "bhtree_worker_sphsum.cpp"
typedef sphlatch::CubicSpline3D                  krnlT;

typedef sphlatch::densSum<partT, krnlT>          densT;
typedef sphlatch::SPHsumWorker<densT, partT>     densSumT;

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>              costT;

#ifdef SPHLATCH_ANEOS
 #include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT>                   eosT;
#else
 #include "eos_idealgas.cpp"
typedef sphlatch::IdealGas<partT>                eosT;
#endif

#include "clump_finder.cpp"

// particles are global
partSetT parts;

void derive()
{
   treeT& Tree(treeT::instance());
   logT&  Logger(logT::instance());

   const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   Logger << "created tree";

   Tree.update(0.8, 1.2);

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   Logger.stream << "Tree.update() -> " << noCZbottomLoc << " CZ cells";
   Logger.flushStream();

   densSumT densWorker(&Tree);
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);
   Logger << "Tree.densWorker()";

   eosT& EOS(eosT::instance());
   for (size_t i = 0; i < nop; i++)
   {
#ifdef SPHLATCH_TIMEDEP_ENERGY
      if ( parts[i].u < uMin )
        parts[i].u = uMin;
#endif
      EOS(parts[i]);
   }
   Logger << "pressure";

   
   costT costWorker(&Tree);
#pragma omp parallel for firstprivate(costWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      costWorker(CZbottomLoc[i]);
   Logger << "Tree.costWorker()";

   Tree.clear();
   Logger << "Tree.clear()";
}

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
   derive();

   getClumps(parts);
   
   parts.doublePrecOut();
   parts.saveHDF5("out.h5part");
   Logger << "stored dump ";

   return(0);
}

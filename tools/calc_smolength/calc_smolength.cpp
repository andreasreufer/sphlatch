#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5

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

class particle :
   public sphlatch::treePart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart
{
public:
   
  ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(h, "h"));
      return(vars);
   }
};

typedef particle   partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#include "bhtree_worker_neighfunc.cpp"
typedef sphlatch::SmoLenFindWorker<partT> smolT;

/*#include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;*/


// particles are global
partSetT parts;

void derive()
{
}

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not ((argc == 3) || (argc == 4)))
   {
      std::cerr <<
      "usage: nbody_onerun_XXXXXX <inputdump> <outputdump> (<numthreads>)\n";
      return(1);
   }


   std::string inFilename = argv[1];
   std::string outFilename = argv[2];

   if (argc == 4)
   {
      std::istringstream threadStr(argv[3]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }

   ///
   /// log program compilation time
   ///
   std::cout  << "executable compiled from " << __FILE__
                 << " on " << __DATE__
                 << " at " << __TIME__ << "\n\n";

   // load the particles
   parts.loadHDF5(inFilename);

   // create tree
   treeT& Tree(treeT::instance());

   const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }

   Tree.update(0.8, 1.2);

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   // estimate h
   smolT       smolWorker(&Tree);
 
#pragma omp parallel for firstprivate(smolWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
     smolWorker(CZbottomLoc[i]);

   Tree.clear();

   parts.doublePrecOut();
   parts.saveHDF5(outFilename);

   return(0);
}

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
typedef sphlatch::BHTree    treeT;

#define SPHLATCH_LOGGER
#include "logger.cpp"
typedef sphlatch::Logger     logT;

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
   fType   accrelerr;
   vect3dT accbf;

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

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
      vars.push_back(storeVar(acc, "acc"));
      vars.push_back(storeVar(accbf, "accbf"));
      vars.push_back(storeVar(accrelerr, "accrelerr"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
#ifdef SPHLATCH_GRAVITY_EPSSMOOTHING
      vars.push_back(storeVar(eps, "eps"));
#endif
      return(vars);
   }
};

typedef particle                               partT;

#include "particle_set.cpp"
typedef sphlatch::ParticleSet<partT>           partSetT;

#ifdef SPHLATCH_GRAVITY
 #include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;
#endif

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

   Tree.update(0.8, 1.2);
   Logger << "created tree";

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   for (size_t i = 0; i < nop; i++)
      parts[i].acc = 0., 0., 0.;

#ifdef SPHLATCH_GRAVITY
   const fType G = parts.attributes["gravconst"];
   gravT       gravWorker(&Tree, G);
 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcAcc(CZbottomLoc[i]);
#endif
   Tree.clear();

   Logger << "Tree gravity";

   std::cout << nop << "\n";
   for (size_t i = 0; i < nop; i++)
   {
      vect3dT accbf;
      accbf = 0., 0., 0.;
      for (size_t j = 0; j < nop; j++)
      {
         if (i != j) 
         {
            const vect3dT rv = parts[i].pos - parts[j].pos;
            const fType   rr = dot(rv, rv);
            const fType   r  = sqrt(rr);
            const fType   mj = parts[j].m;

#ifdef SPHLATCH_GRAVITY_SPLINESMOOTHING
            const fType hj   = parts[j].h;
            const fType mOr3 = mj * gravWorker.splineOSmoR3(r, hj);
#elif SPHLATCH_GRAVITY_EPSSMOOTHING
            const fType epsj = parts[j].eps;
            const fType re   = r + epsj;
            const fType mOr3 = mj / (re * re * re);
#else
            const fType mOr3 = mj / (rr * r);
#endif
            accbf -= mOr3 * rv;
         }
      }
      parts[i].accbf = accbf * G;

      const vect3dT dacc  = parts[i].acc - parts[i].accbf;
      const fType   daccs = sqrt(dot(dacc, dacc));
      const fType   accs  = sqrt(dot(parts[i].acc, parts[i].acc));

      parts[i].accrelerr = daccs / accs;
   }
   Logger << "BF   gravity";
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


   std::string inFilename  = argv[1];
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
   std::cout << "executable compiled from " << __FILE__
             << " on " << __DATE__
             << " at " << __TIME__ << "\n\n";

   // load the particles
   parts.loadHDF5(inFilename);

   derive();

   parts.doublePrecOut();
   parts.saveHDF5(outFilename);

   return(0);
}

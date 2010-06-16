#include <iostream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

//#include <omp.h>
//#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
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
#include "integrator_predcorr.cpp"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::energyPart,
   public sphlatch::IOPart
#ifdef SPHLATCH_TIMEDEP_SMOOTHING
   , public sphlatch::varHPart
#endif
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
      vars.push_back(storeVar(acc, "acc"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(rho, "rho"));
      vars.push_back(storeVar(noneigh, "noneigh"));

      vars.push_back(storeVar(p, "p"));
      vars.push_back(storeVar(u, "u"));
      vars.push_back(storeVar(cs, "cs"));
#ifdef SPHLATCH_ANEOS
      vars.push_back(storeVar(mat, "mat"));
      vars.push_back(storeVar(T, "T"));
      vars.push_back(storeVar(phase, "phase"));
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

#ifdef SPHLATCH_GRAVITY
 #include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;
#endif

#include "sph_algorithms.cpp"
#include "sph_kernels.cpp"
#include "bhtree_worker_sphsum.cpp"

typedef sphlatch::CubicSpline3D                  krnlT;

typedef sphlatch::densSum<partT, krnlT>          densT;
typedef sphlatch::SPHsumWorker<densT, partT>     densSumT;

typedef sphlatch::SPHsumWorker<densT, partT>     densSum2T;

//typedef sphlatch::accPowSum<partT, krnlT>        accPowT;
//typedef sphlatch::SPHsumWorker<accPowT, partT>   accPowSumT;

#ifdef SPHLATCH_ANEOS
#include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT>                   eosT;
#endif

// particles are global
partSetT parts;

void derive()
{
   treeT& Tree(treeT::instance());
   logT&  Logger(logT::instance());

   Tree.update(0.8, 1.2);
   Logger << "Tree.update()";

   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   const size_t nop = parts.getNop();
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].acc  = 0., 0., 0.;
      parts[i].dudt = 0.;
      parts[i].rho = 0.;
   }
/*
#ifdef SPHLATCH_GRAVITY
   const fType G = parts.attributes["gravconst"];
   gravT gravWorker(&Tree, G);
 #pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcAcc(CZbottomLoc[i]);
   Logger << "Tree.calcAcc()";
#endif
  */ 
   densSumT densWorker(&Tree);
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);
   Logger << "Tree.densWorker()";
   
   /*eosT& EOS(eosT::instance());
   for (size_t i = 0; i < nop; i++)
   {
      EOS(parts[i]);
   }
   Logger << "pressure";*/
/*
   for (size_t i = 0; i < nop; i++)
   {
      parts[i].p = 1.e5;
   }

   accPowSumT accPowWorker(&Tree);
#pragma omp parallel for firstprivate(accPowWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      accPowWorker(CZbottomLoc[i]);
   Logger << "Tree.accPowWorker()";
   */
}

int main(int argc, char* argv[])
{
   logT& Logger(logT::instance());

   if ( argc != 2 )
     return(1);

   std::string inFilename = argv[1];

   // load the particles
   parts.loadHDF5(inFilename);
   Logger << "loaded particles";

   treeT& Tree(treeT::instance());

   const size_t nop       = parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(parts.getBox() * 1.5);

   for (size_t i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   Logger << "created tree";

   // start the loop
   for (size_t i = 0; i < 1; i++)
   {
      derive();
   }

   return(0);
}

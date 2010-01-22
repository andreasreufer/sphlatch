#include <iostream>
#include <sstream>
#include <vector>

//#define SPHLATCH_SINGLEPREC

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#ifdef SPHLATCH_TIMEDEP_ENERGY
 #ifndef SPHLATCH_VELDIV
  #define SPHLATCH_VELDIV
 #endif
#endif

#ifdef SPHLATCH_TIMEDEP_SMOOTHING
 #ifndef SPHLATCH_VELDIV
  #define SPHLATCH_VELDIV
 #endif
 #ifndef SPHLATCH_NONEIGH
  #define SPHLATCH_NONEIGH
 #endif
#endif

#ifdef SPHLATCH_INTEGRATERHO
 #ifndef SPHLATCH_VELDIV
  #define SPHLATCH_VELDIV
 #endif
#endif

#ifdef SPHLATCH_ANEOS
 #ifndef SPHLATCH_NONEGPRESS
  #define SPHLATCH_NONEGPRESS
 #endif
#endif


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

      // minimal
      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(u, "u"));

      if (not fName.empty())
         vars.push_back(storeVar(f, fName));

      if (not vName.empty())
         vars.push_back(storeVar(v, vName));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      return(vars);
   }

   static std::string fName, vName;

public:
   fType   f, ftar;
   vect3dT v, vtar;

   void setLoadF(std::string _str)
   {
      fName = _str;
   }

   void setLoadV(std::string _str)
   {
      vName = _str;
   }
};

std::string particle::fName;
std::string particle::vName;

template<typename _partT, typename _krnlT>
struct interpolSum
{
   _krnlT  K;

   void    preSum(_partT* const _i)
   {}

   void operator()(_partT* const _i,
                   _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      const fType r  = sqrt(_rr);
      const fType hj = _j->h;

      const fType rhoi = _i->rho;
      const fType mi   = _i->m;

      const fType k = (K.value(r, hj)) * (mi / rhoi);

      _j->ftar += k * _i->f;
      _j->vtar += k * _i->v;
   }

   void postSum(_partT* const _i)
   {}
};

typedef particle                               partT;

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

typedef sphlatch::CubicSpline3D                    krnlT;

typedef sphlatch::densSum<partT, krnlT>          densT;
typedef sphlatch::SPHsumWorker<densT, partT>     densSumT;

typedef interpolSum<partT, krnlT>        interpolSumT;
typedef sphlatch::SPHsumWorker<interpolSumT, partT>   interpolWorkT;

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>                costT;

#ifdef SPHLATCH_ANEOS
 #include "eos_aneos.cpp"
typedef sphlatch::ANEOS<partT>                     eosT;
#else
 #include "eos_idealgas.cpp"
typedef sphlatch::IdealGas<partT>                  eosT;
#endif


// particles are global
partSetT parts;
partSetT gridParts;



void derive()
{
   treeT& Tree(treeT::instance());

   const int nop       = parts.getNop();
   const size_t nogp       = gridParts.getNop();
   
   const fType  costppart = 1. / ( nop + nogp);

   Tree.setExtent(parts.getBox() * 1.1);

   for (int i = 0; i < nop; i++)
   {
      parts[i].cost = costppart;
      Tree.insertPart(parts[i]);
   }
   
   for (size_t i = 0; i < nogp; i++)
   {
      gridParts[i].cost = costppart;
      gridParts[i].m = 0.;
      gridParts[i].rho = 1.;
      Tree.insertPart(gridParts[i]);
   }
   std::cout << "created tree";

   Tree.update(0.8, 1.2);

   densSumT densWorker(&Tree);
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < nop; i++)
      densWorker(&(parts[i]));
   std::cout << "Tree.densWorker()";

   interpolWorkT interpolWorker(&Tree);
#pragma omp parallel for firstprivate(interpolWorker)
   for (int i = 0; i < nop; i++)
      interpolWorker(&(parts[i]));
   std::cout << "Tree.interpolWorker()";

   Tree.clear();
   std::cout << "Tree.clear()";
}

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not ((argc == 3) || (argc == 4)))
   {
      std::cerr <<
      "usage: sph2grid <inputdump> <outputFile> (<numthreads>)\n";
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

   // load the particles
   parts.loadHDF5(inFilename);
   std::cout << "loaded particles";

   fType& time(parts.attributes["time"]);

   gridParts.step = parts.step;
   gridParts.attributes = parts.attributes;

   // first bootstrapping step
   derive();

   gridParts.doublePrecOut();
   gridParts.saveHDF5(outFilename);

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}


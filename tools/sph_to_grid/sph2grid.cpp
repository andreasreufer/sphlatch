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
typedef sphlatch::fmatrT    fmatrT;

const fType finf = sphlatch::fTypeInf;

#include "parse_po_vectors.cpp"

#include "bhtree.cpp"
typedef sphlatch::BHTree   treeT;

#include "hdf5_io.cpp"
typedef sphlatch::HDF5File H5FT;

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

#ifdef SPHLATCH_INTERPOLATE_SCALAR
      if (not fName.empty())
         vars.push_back(storeVar(f, fName));
#endif

#ifdef SPHLATCH_INTERPOLATE_VECTOR
      if (not vName.empty())
         vars.push_back(storeVar(v, vName));
#endif

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;
      return(vars);
   }


public:
#ifdef SPHLATCH_INTERPOLATE_SCALAR
   static std::string fName;
   
   fType   f, ftar;
   
   void setLoadF(std::string _str)
   {
      fName = _str;
   }
#endif


#ifdef SPHLATCH_INTERPOLATE_VECTOR
   static std::string vName;
   
   vect3dT v, vtar;

   void setLoadV(std::string _str)
   {
      vName = _str;
   }
#endif
};

#ifdef SPHLATCH_INTERPOLATE_SCALAR
std::string particle::fName;
#endif
#ifdef SPHLATCH_INTERPOLATE_VECTOR
std::string particle::vName;
#endif

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
      const fType hi = _i->h;

#ifdef SPHLATCH_MISCIBLE
      const fType deltai = _i->delta;
      const fType k = (K.value(r, hi)) * (1./deltai);
#else
      const fType rhoi = _i->rho;
      const fType mi   = _i->m;
      const fType k = (K.value(r, hi)) * (mi / rhoi);
#endif

#ifdef SPHLATCH_INTERPOLATE_SCALAR
      _j->ftar += k * _i->f;
#endif
      
#ifdef SPHLATCH_INTERPOLATE_VECTOR
      _j->vtar += k * _i->v;
#endif
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

   Tree.setExtent(( parts.getBox() + gridParts.getBox() ) * 1.1);

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
   std::cout << "created tree\n";

   Tree.update(0.8, 1.2);

   densSumT densWorker(&Tree);
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < nop; i++)
      densWorker(&(parts[i]));
   std::cout << "Tree.densWorker()\n";

   interpolWorkT interpolWorker(&Tree);
#pragma omp parallel for firstprivate(interpolWorker)
   for (int i = 0; i < nop; i++)
      interpolWorker(&(parts[i]));
   std::cout << "Tree.interpolWorker()\n";

   Tree.clear();
   std::cout << "Tree.clear()\n";
}

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif

   if (not ((argc == 9) || (argc == 10)))
   {
      std::cerr <<
      "usage: sph2grid_X <inputdump> <outputFile> <var> <nx> <ny> <r0> <ri> <rj> (<numthreads>)\n";
      return(1);
   }

   std::string inFilename = argv[1];
   std::string outFilename = argv[2];
   std::string varname     = argv[3];

   std::istringstream niStr(argv[4]);
   size_t ni;
   niStr >> ni;
   
   std::istringstream njStr(argv[5]);
   size_t nj;
   njStr >> nj;

   
   const vect3dT r0 = vectOptParse(argv[6]);
   const vect3dT ri = vectOptParse(argv[7]);
   const vect3dT rj = vectOptParse(argv[8]);
      

   if (argc == 10)
   {
      std::istringstream threadStr(argv[9]);
      int numThreads;
      threadStr >> numThreads;
      omp_set_num_threads(numThreads);
   }
  
   parts.resize(1);
   
#ifdef SPHLATCH_INTERPOLATE_SCALAR
   parts[0].setLoadF(varname);
#endif
#ifdef SPHLATCH_INTERPOLATE_VECTOR
   parts[0].setLoadV(varname);
#endif
   
   // load the particles
   parts.loadHDF5(inFilename);
   std::cout << "loaded particles\n";

   gridParts.step = parts.step;
   gridParts.attributes = parts.attributes;

   const size_t nogp = ni*nj;
   gridParts.resize(nogp);

   for (size_t i = 0; i < ni; i++)
     for (size_t j = 0; j < nj; j++)
     {
       const size_t k = j*ni + i;
       const fType fi = i;
       const fType fj = j;
       gridParts[k].pos = r0 + ri*fi + rj*fj;
       
       gridParts[k].id = -1;
       gridParts[k].m = 0.;
       gridParts[k].h = 0.;
       gridParts[k].rho = 1.;
     }

  
   // interpolate
   derive();

   fmatrT ripos(ni,3), rjpos(nj,3);

   for (size_t i = 0; i < ni; i++)
   {
     ripos(i, 0) = r0[0] + i*ri[0];
     ripos(i, 1) = r0[1] + i*ri[1];
     ripos(i, 2) = r0[2] + i*ri[2];
   }
   
   for (size_t j = 0; j < nj; j++)
   {
     ripos(j, 0) = r0[0] + j*rj[0];
     ripos(j, 1) = r0[1] + j*rj[1];
     ripos(j, 2) = r0[2] + j*rj[2];
   }
   
   
#ifdef SPHLATCH_INTERPOLATE_SCALAR
   fmatrT f(ni, nj);
#endif

#ifdef SPHLATCH_INTERPOLATE_VECTOR
   fmatrT v0(ni,nj),  v1(ni,nj), v2(ni,nj);
#endif

   for (size_t i = 0; i < ni; i++)
     for (size_t j = 0; j < nj; j++)
     {
       const size_t k = j*ni + i;

#ifdef SPHLATCH_INTERPOLATE_SCALAR
       f(i,j)= gridParts[k].ftar;
#endif

#ifdef SPHLATCH_INTERPOLATE_VECTOR
       v0(i,j) =  gridParts[k].vtar(0);
       v1(i,j) =  gridParts[k].vtar(1);
       v2(i,j) =  gridParts[k].vtar(2);
#endif
     }
   
   H5FT gridFile(outFilename);
   
   gridFile.savePrimitive("ripos",ripos);
  
#ifdef SPHLATCH_INTERPOLATE_SCALAR
   gridFile.savePrimitive(varname,f);
#endif
   
#ifdef SPHLATCH_INTERPOLATE_VECTOR
   gridFile.savePrimitive(varname + "_x",v0);
   gridFile.savePrimitive(varname + "_y",v1);
   gridFile.savePrimitive(varname + "_z",v2);
#endif

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}


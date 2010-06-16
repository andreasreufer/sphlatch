#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::vect3dT   vect3dT;
typedef sphlatch::box3dT    box3dT;

#include "bhtree.cpp"
typedef sphlatch::BHTree    treeT;

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
{
public:

   sphlatch::PredictorCorrectorO2<vect3dT> posInt;
   sphlatch::PredictorCorrectorO1<fType>   energyInt;

   void bootstrap()
   {
      posInt.bootstrap(pos, vel, acc);
      energyInt.bootstrap(u, dudt);
   }

   void predict(const fType _dt)
   {
      posInt.predict(pos, vel, acc, _dt);
      energyInt.predict(u, dudt, _dt);
   }

   void correct(const fType _dt)
   {
      posInt.correct(pos, vel, acc, _dt);
      energyInt.correct(u, dudt, _dt);
   }

   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(h, "h"));
      vars.push_back(storeVar(id, "id"));

      vars.push_back(storeVar(u, "u"));

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

      vars.push_back(storeVar(cost, "cost"));

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

#include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

#include "sph_algorithms.cpp"
#include "sph_kernels.cpp"
#include "bhtree_worker_sphsum.cpp"

typedef sphlatch::CubicSpline3D                  krnlT;

typedef sphlatch::densSum<partT, krnlT>          densT;
typedef sphlatch::SPHsumWorker<densT, partT>     densSumT;

typedef sphlatch::accPowSum<partT, krnlT>        accPowT;
typedef sphlatch::SPHsumWorker<accPowT, partT>   accPowSumT;

#include "bhtree_worker_cost.cpp"
typedef sphlatch::CostWorker<partT>              costT;

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif
   treeT& Tree(treeT::instance());

   partSetT partsTree;
   partsTree.loadHDF5("in.h5part");

   const size_t nop    = partsTree.getNop();
   const fType  cppart = 1. / nop;

   std::cout << nop << " " << cppart << "\n";

   Tree.setExtent(partsTree.getBox() * 1.5);

   double start;
   start = omp_get_wtime();
   for (size_t i = 0; i < nop; i++)
   {
      partsTree[i].cost = cppart;
      Tree.insertPart(partsTree[i]);
   }
   std::cout << "partsTree insert " << omp_get_wtime() - start << "s\n";

   partSetT partsBF = partsTree;

   start = omp_get_wtime();
   Tree.update(0.8, 1.2);
   std::cout << "Tree.update()    " << omp_get_wtime() - start << "s\n";

   const fType         G = 1.;
   gravT               gravWorker(&Tree, G);
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   start = omp_get_wtime();
#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcAcc(CZbottomLoc[i]);

   std::cout << "calcAcc()    " << omp_get_wtime() - start << "s\n";

   densSumT densWorker(&Tree);
   start = omp_get_wtime();
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);
   std::cout << "densWorker()     " << omp_get_wtime() - start << "s\n";

   accPowSumT accPowWorker(&Tree);
   start = omp_get_wtime();
#pragma omp parallel for firstprivate(accPowWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      accPowWorker(CZbottomLoc[i]);
   std::cout << "accPowWorker()   " << omp_get_wtime() - start << "s\n";

   costT costWorker(&Tree);
   start = omp_get_wtime();
#pragma omp parallel for firstprivate(costWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      costWorker(CZbottomLoc[i]);
   std::cout << "costWorker()     " << omp_get_wtime() - start << "s\n";


   std::cout << "number of CZ cells was " << noCZbottomLoc << "\n";
   partsTree.saveHDF5("out_tree.h5part");
   std::cout << "particle size: " << sizeof(partT) << "\n";
   
   fType totCost = 0.;
   for (size_t i = 0; i < nop; i++)
      totCost += partsTree[i].cost;

   std::cout << "total cost was " << totCost << "\n";

#ifdef BRUTEFORCE
   vect3dT cacc;
   size_t  non;
   start = omp_get_wtime();
   for (size_t i = 0; i < nop; i++)
   {
      cacc = 0, 0, 0;
      non  = 1;

      const fType srad = 2. * partsBF[i].h;

      for (size_t j = 0; j < nop; j++)
      {
         if (i != j)
         {
            const vect3dT rvec = partsBF[i].pos - partsBF[j].pos;
            const fType   r    = sqrt(rvec[0] * rvec[0] +
                                      rvec[1] * rvec[1] +
                                      rvec[2] * rvec[2]);
            const fType mOr3 = partsBF[j].m / (r * r * r);

            cacc -= mOr3 * rvec;

            if (r < srad)
               non++;
         }
      }
      partsBF[i].acc     = cacc;
      partsBF[i].noneigh = non;
   }
   std::cout << "calcAcc() BF " << omp_get_wtime() - start << "s\n";
   partsBF.saveHDF5("out_bf.h5part");

   std::cout << "brute force: " << partsBF[9].acc << "\n";
   std::cout << "tree       : " << partsTree[9].acc << "\n";

   size_t nonDiffSum = 0;
   for (size_t i = 0; i < nop; i++)
      nonDiffSum += abs(partsTree[i].noneigh - partsBF[i].noneigh);

   std::cout << " sum of difference in neighbours " << nonDiffSum << "\n";
#endif

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

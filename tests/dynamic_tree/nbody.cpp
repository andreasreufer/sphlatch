#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5
#define SPHLATCH_NONEIGH

#include "typedefs.h"
typedef sphlatch::fType     fType;
typedef sphlatch::vect3dT   vect3dT;

#include "bhtree.cpp"
typedef sphlatch::BHTree    treeT;

#include "bhtree_particle.h"
#include "sph_fluid_particle.h"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart,
   public sphlatch::IOPart
{
public:
   ioVarLT getLoadVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(id, "id"));

      return(vars);
   }

   ioVarLT getSaveVars()
   {
      ioVarLT vars;

      vars.push_back(storeVar(pos, "pos"));
      vars.push_back(storeVar(vel, "vel"));
      vars.push_back(storeVar(vel, "acc"));
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(id, "id"));
      vars.push_back(storeVar(noneigh, "noneigh"));

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
typedef sphlatch::densSum<partT>               densT;
#include "bhtree_worker_sphsum.cpp"
typedef sphlatch::SPHsumWorker<densT, partT>   densSumT;


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

   double start;

   start = omp_get_wtime();
   for (size_t i = 0; i < nop; i++)
   {
      partsTree[i].cost = cppart;
      partsTree[i].h    = 0.1;
      Tree.insertPart(partsTree[i]);
   }
   std::cout << "partsTree insert " << omp_get_wtime() - start << "s\n";

   partSetT partsBF = partsTree;

   start = omp_get_wtime();
   Tree.update(0.8, 1.2);
   std::cout << "Tree.update()    " << omp_get_wtime() - start << "s\n";

   gravT gravWorker(&Tree);
   std::cout << "gravity ...\n";
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   start = omp_get_wtime();
#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcGravity(CZbottomLoc[i]);

   std::cout << "calcGravity()    " << omp_get_wtime() - start << "s\n";

   densSumT densWorker(&Tree);
   start = omp_get_wtime();
#pragma omp parallel for firstprivate(densWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      densWorker(CZbottomLoc[i]);

   std::cout << "densWorker()     " << omp_get_wtime() - start << "s\n";

   partsTree.saveHDF5("out_tree.h5part");

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
   std::cout << "calcGravity() BF " << omp_get_wtime() - start << "s\n";
   partsBF.saveHDF5("out_bf.h5part");

   std::cout << "brute force: " << partsBF[9].acc << "\n";
   std::cout << "tree       : " << partsTree[9].acc << "\n";

   size_t nonDiffSum = 0;
   for (size_t i = 0; i < nop; i++)
      nonDiffSum += abs(partsTree[i].noneigh - partsBF[i].noneigh);

   std::cout << " sum of difference in neighbours " << nonDiffSum << "\n";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

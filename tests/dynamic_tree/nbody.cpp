#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP
#define SPHLATCH_HDF5

#include "typedefs.h"
typedef sphlatch::fType    fType;
typedef sphlatch::vect3dT  vect3dT;

#include "bhtree.cpp"
typedef sphlatch::BHTree   treeT;

#include "bhtree_particle.h"
#include "io_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
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
      vars.push_back(storeVar(m, "m"));
      vars.push_back(storeVar(id, "id"));

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

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif
   treeT& Tree(treeT::instance());

   partSetT particles;
   particles.loadHDF5("in.h5part");
   const size_t nop = particles.getNop();
   const fType cppart = 1. / nop;
   std::cout << nop << " " << cppart << "\n";
   
   double start;

   start = omp_get_wtime();
   for (size_t i = 0; i < nop; i++)
   {
     particles[i].cost = cppart;
     Tree.insertPart(particles[i]);
   }
   std::cout << "particles insert " << omp_get_wtime() - start << "s\n";

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

   particles.saveHDF5("out_tree.h5part");

   partSetT partsBF;
   partsBF.loadHDF5("in.h5part");

   vect3dT cacc;
   start = omp_get_wtime();
   for (size_t i = 0; i < nop; i++)
   {
      cacc = 0, 0, 0;
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
         }
      }
      partsBF[i].acc = cacc;
   }
   std::cout << "calcGravity() BF " << omp_get_wtime() - start << "s\n";
   partsBF.saveHDF5("out_bf.h5part");

   std::cout << "brute force: " << partsBF[9].acc << "\n";
   std::cout << "tree       : " << particles[9].acc << "\n";
  
   /*particles[9].acc = 0,0,0;
   gravWorker.calcGravPartAlt(particles[9].treeNode);
   std::cout << "\n";
   std::cout << particles[9].acc << "\n\n";*/

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

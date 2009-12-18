#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP

#include "typedefs.h"
typedef sphlatch::fType      fType;

#include "bhtree.cpp"
typedef sphlatch::BHTree     treeT;

typedef sphlatch::pnodT      pnodT;
typedef sphlatch::pnodPtrT   pnodPtrT;

typedef sphlatch::nodeT      nodeT;
typedef sphlatch::nodePtrT   nodePtrT;

typedef sphlatch::gcllT      gcllT;
typedef sphlatch::gcllPtrT   gcllPtrT;

typedef sphlatch::czllT      czllT;
typedef sphlatch::czllPtrT   czllPtrT;


#include "bhtree_particle.h"
#include "sph_fluid_particle.h"

class particle :
   public sphlatch::treePart,
   public sphlatch::movingPart,
   public sphlatch::SPHfluidPart { };

class ghost :
   public sphlatch::treeGhost,
   public sphlatch::movingGhost,
   public sphlatch::SPHfluidGhost { };

typedef particle                               partT;
typedef ghost                                  ghstT;

#include "bhtree_treedump.cpp"
typedef sphlatch::BHTreeDump                   dumpT;

#include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

//#include "log_manager.h"
//typedef sphlatch::LogManager                   logT;

struct densFunc
{
   void operator()(partT* _i, const partT* _j)
   {
      _i->pos = _j->pos;
   }
};

#include "bhtree_worker_sphsum.cpp"
typedef sphlatch::SPHsumWorker<densFunc, partT>   densSumT;

#include "bhtree_worker_neighfunc.cpp"
typedef sphlatch::NeighWorker<densFunc, partT> neighFuncT;


#include "bhtree_worker.h"
class BHTreeTester : public sphlatch::BHTreeWorker {
public:
   BHTreeTester(treePtrT _treePtr)
      : BHTreeWorker(_treePtr) { }
   ~BHTreeTester() { }

   void build()
   { }

   void test()
   { }

   void dispRoot()
   {
      std::cout << static_cast<czllPtrT>(rootPtr)->noParts << "\n";
      std::cout << static_cast<czllPtrT>(rootPtr)->relCost << "\n";
   }
};

int main(int argc, char* argv[])
{
#ifdef SPHLATCH_MPI
   MPI::Init(argc, argv);
#endif
   treeT& Tree(treeT::instance());
   dumpT  dumper(&Tree);

   BHTreeTester testWorker(&Tree);

   const fType costMin = 1.0e4, costMax = 1.5e4;
   const size_t       noParts = 1000000;
   //const size_t noParts = 300000;

   //const fType costMin = 10., costMax = 15.;
   //const size_t noParts = 100;

   std::vector<partT> particles(noParts);
   for (size_t i = 0; i < noParts; i++)
   {
      particles[i].pos[0] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[1] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[2] = static_cast<fType>(rand()) / RAND_MAX;

      particles[i].m = 1.;

      particles[i].id   = i;
      particles[i].cost = 1.;
   }

   double start;
   start = omp_get_wtime();
   for (size_t i = 0; i < noParts; i++)
   {
      Tree.insertPart(particles[i]);
   }
   std::cout << "particles insert " << omp_get_wtime() - start << "s\n";

   std::cout << "Tree.update() .........\n";
   start = omp_get_wtime();
   Tree.update(costMin, costMax);
   std::cout << "Tree.update()    " << omp_get_wtime() - start << "s\n";

   //dumper.ptrDump("postCZdump.txt");
   //dumper.dotDump("postCZdump.dot");

   /*testWorker.dispRoot();
      testWorker.build();
      testWorker.test();*/

   /*for (size_t i = 0; i < noParts; i++)
   {
      if (particles[i].pos[2] < 0.5)
         particles[i].pos[2] += 0.25;

      particles[i].cost = 1.;

   }*/

   std::cout << "Tree.update() .........\n";
   start = omp_get_wtime();
   Tree.update(costMin, costMax);
   std::cout << "Tree.update()    " << omp_get_wtime() - start << "s\n";
   
   //dumper.dotDump("post_move.dot");
   //dumper.ptrDump("post_move.txt");

   const fType G = 1.;
   gravT gravWorker(&Tree, G);

   std::cout << "gravity ...\n";
   treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   const int           noCZbottomLoc = CZbottomLoc.size();

   start = omp_get_wtime();
#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
      gravWorker.calcGravity(CZbottomLoc[i]);

   std::cout << "gravity finished " << omp_get_wtime() - start << "s\n";

   //for (size_t i = 0; i < noParts; i++)
   //std::cout << particles[i].acc << "\n";

#ifdef SPHLATCH_MPI
   MPI::Finalize();
#endif
   return(0);
}

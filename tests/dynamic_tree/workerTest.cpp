#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP

#include "typedefs.h"
typedef sphlatch::fType        fType;

#include "bhtree.cpp"
typedef sphlatch::BHTree       treeT;

typedef sphlatch::pnodT        pnodT;
typedef sphlatch::pnodPtrT     pnodPtrT;

typedef sphlatch::nodeT        nodeT;
typedef sphlatch::nodePtrT     nodePtrT;

typedef sphlatch::gcllT        gcllT;
typedef sphlatch::gcllPtrT     gcllPtrT;

typedef sphlatch::czllT        czllT;
typedef sphlatch::czllPtrT     czllPtrT;

#include "bhtree_treedump.cpp"
typedef sphlatch::BHTreeDump   dumpT;

#include "bhtree_particle.h"
typedef sphlatch::treeGhost    partT;

struct densFunc
{
   void operator()(partT* _i, const partT* _j)
   {
      _i->pos = _j->pos;
   }
};

#include "bhtree_worker_sphsum.cpp"
typedef sphlatch::SPHsumWorker<densFunc>   densSumT;

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
      std::cout << static_cast<czllPtrT>(rootPtr)->absCost << "\n";
   }
};

int main(int argc, char* argv[])
{
   sleep(1);
   MPI::Init(argc, argv);

   treeT&       Tree(treeT::instance());
   dumpT        dumper(&Tree);
   BHTreeTester testWorker(&Tree);

   const size_t       noParts = 100;
   std::vector<partT> particles(noParts);
   for (size_t i = 0; i < noParts; i++)
   {
      particles[i].pos[0] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[1] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[2] = static_cast<fType>(rand()) / RAND_MAX;

      particles[i].id   = i;
      particles[i].cost = 1.;
   }

   Tree.insertParts(particles);
   Tree.update();

   //dumper.ptrDump("postCZdump.txt");
   //dumper.dotDump("postCZdump.dot");

   testWorker.dispRoot();
   testWorker.build();
   testWorker.test();

   for (size_t i = 0; i < noParts; i++)
   {
      if (particles[i].pos[0] < 0.5)
         particles[i].pos[0] += 0.5;

      particles[i].cost = 1.;
   }
   //Tree.update();

   MPI::Finalize();
   return(0);
}

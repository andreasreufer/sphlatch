#include <iostream>
#include <vector>

#define SPHLATCH_OPENMP
#include <omp.h>

#include "typedefs.h"

//#include "particle_manager.h"
//typedef sphlatch::ParticleManager          partMgrT;

#include "bhtree_dynamic.h"
typedef sphlatch::BHTree                   treeT;
//typedef sphlatch::BHTree::BHTreeConfig     treeConfigT;

typedef sphlatch::particleNode             pnodT;

typedef sphlatch::genericNode              nodeT;
typedef sphlatch::genericNode*             nodePtrT;

typedef sphlatch::quadrupoleCellNode       cellT;
typedef sphlatch::quadrupoleCellNode*      cellPtrT;

typedef sphlatch::costzoneCellNode         czllT;

#include "bhtree_part_insertmover.h"
typedef sphlatch::BHTreePartsInsertMover   inserterT;

#include "bhtree_treedump.h"
typedef sphlatch::BHTreeDump               dumpT;

#include "particle.h"
typedef sphlatch::treeParticle             partT;

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
};

int main(int argc, char* argv[])
{
   sleep(1);
   //MPI::Init(argc, argv);

   treeT& Tree(treeT::instance());

   inserterT    inserter(&Tree);
   dumpT        dumper(&Tree);
   BHTreeTester testWorker(&Tree);
   std::cout << "tree workers instantiated\n";


   std::cout << "push down particle ...\n";

   partT part1, part2;

   part1.pos = 0.,0.,0.;
   part1.id = 0;

   part2.pos = 1.,1.,1.;
   part2.id = 1;
   
   std::cout << "insert particle ...\n";
   inserter.insert(part1);
   inserter.insert(part2);
   std::cout << " done!\n";
   dumper.dotDump("test0.dot");
   //dumper.ptrDump();

   inserter.pushDown(part1);
   //dumper.ptrDump();
   inserter.pushDown(part2);
   //dumper.ptrDump();
   std::cout << " done!\n";
   dumper.dotDump("test1.dot");

   //dumper.ptrDump();

   /*testWorker.build();
      testWorker.test();*/

   //MPI::Finalize();
   return(0);
}

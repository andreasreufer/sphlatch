#include <iostream>
#include <vector>

#define SPHLATCH_OPENMP
#include <omp.h>

#include "typedefs.h"
typedef sphlatch::matrixRefType            matrixRefT;

//#include "particle_manager.h"
//typedef sphlatch::ParticleManager          partMgrT;

#include "bhtree_dynamic.h"
typedef sphlatch::BHTree                   treeT;
//typedef sphlatch::BHTree::BHTreeConfig     treeConfigT;

typedef sphlatch::particleNode             partT;

typedef sphlatch::genericNode              nodeT;
typedef sphlatch::genericNode*             nodePtrT;

typedef sphlatch::quadrupoleCellNode       cellT;
typedef sphlatch::quadrupoleCellNode*      cellPtrT;

typedef sphlatch::costzoneCellNode         czllT;

#include "bhtree_part_insertmover.h"
typedef sphlatch::BHTreePartsInsertMover   inserterT;

#include "bhtree_treedump.h"
typedef sphlatch::BHTreeDump               dumpT;

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

   partMgrT&  PartManager(partMgrT::instance());
   matrixRefT pos(PartManager.pos);


   treeT&      Tree(treeT::instance());

   inserterT    inserter(&Tree);
   dumpT        dumper(&Tree);
   BHTreeTester testWorker(&Tree);
   std::cout << "tree workers instantiated\n";

   std::cout << "insert particle ...\n";
   inserter.insert(0);
   inserter.insert(1);
   std::cout << " done!\n";
   dumper.dotDump("test0.dot");
   //dumper.ptrDump();
   
   std::cout << "push down particle ...\n";
   inserter.pushDown(0);
   //dumper.ptrDump();
   inserter.pushDown(1);
   //dumper.ptrDump();
   std::cout << " done!\n";
   dumper.dotDump("test1.dot");
   //dumper.ptrDump();

   /*testWorker.build();
   testWorker.test();*/

   //MPI::Finalize();
   return(0);
}

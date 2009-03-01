#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP

#include "typedefs.h"
typedef sphlatch::fType                    fType;

#include "bhtree_dynamic.h"
typedef sphlatch::BHTree                   treeT;

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

//#include "bhtree_cz_builder.h"
//typedef sphlatch::BHTreeCZBuilder          czbldT;

#include "particle.h"
typedef sphlatch::treeParticle   partT;

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
   //czbldT       CZbuilder(&Tree);
   std::cout << "tree workers instantiated\n";


   std::cout << "push down particle ...\n";

   const size_t noParts = 100;

   std::vector<partT> particles(noParts);

   for (size_t i = 0; i < noParts; i++)
   {
      particles[i].pos[0] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[1] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[2] = static_cast<fType>(rand()) / RAND_MAX;

      particles[i].id = i;

      inserter.insert(particles[i]);
   }


  //dumper.dotDump("test0.dot");
   for (size_t i = 0; i < noParts; i++)
   {
      inserter.pushDown(particles[i]);

      /*if ( i == 7 )
        dumper.dotDump("test0.dot");
      
      if ( i == 8 )
        dumper.dotDump("test1.dot");*/

   }
  dumper.ptrDump();
  dumper.dotDump("test2.dot");
  //dumper.dotDump("test1.dot");

   /*testWorker.build();
      testWorker.test();*/

   //MPI::Finalize();
   return(0);
}

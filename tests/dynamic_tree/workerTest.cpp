#include <iostream>
#include <vector>

#include <omp.h>
#define SPHLATCH_OPENMP

#include "typedefs.h"
typedef sphlatch::fType                    fType;

#include "bhtree.cpp"
typedef sphlatch::BHTree                   treeT;

typedef sphlatch::pnodT                    pnodT;
typedef sphlatch::pnodPtrT                 pnodPtrT;

typedef sphlatch::nodeT                    nodeT;
typedef sphlatch::nodePtrT                 nodePtrT;

typedef sphlatch::gcllT                    gcllT;
typedef sphlatch::gcllPtrT                 gcllPtrT;

typedef sphlatch::czllT                    czllT;
typedef sphlatch::czllPtrT                 czllPtrT;

#include "bhtree_part_insertmover.cpp"
typedef sphlatch::BHTreePartsInsertMover   inserterT;

//#include "bhtree_treedump.h"
//typedef sphlatch::BHTreeDump               dumpT;

//#include "bhtree_cz_builder.h"
//typedef sphlatch::BHTreeCZBuilder          czbldT;

#include "bhtree_housekeeper.cpp"

#include "bhtree_particle.h"
typedef sphlatch::treeGhost   partT;

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
   //MPI::Init(argc, argv);

   treeT& Tree(treeT::instance());

   inserterT inserter(&Tree);
   //dumpT        dumper(&Tree);
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
      particles[i].cost = 1.;

      std::cout << "try to insert :" << particles[i].pos << "\n";
      inserter.insert(particles[i]);
   }


   //dumper.dotDump("test0.dot");
   for (size_t i = 0; i < noParts; i++)
   {
      std::cout << "push down particle " << i << "\n";
      inserter.pushDown(particles[i]);

      /*if ( i == 7 )
         dumper.dotDump("test0.dot");

         if ( i == 8 )
         dumper.dotDump("test1.dot");*/
   }

   testWorker.dispRoot();

   //dumper.ptrDump();
   //dumper.dotDump("test2.dot");
   //dumper.dotDump("test1.dot");

   /*testWorker.build();
      testWorker.test();*/

   //MPI::Finalize();
   return(0);
}

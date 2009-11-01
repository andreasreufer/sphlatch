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
  public sphlatch::SPHfluidPart {};

class ghost : 
  public sphlatch::treeGhost,  
  public sphlatch::movingGhost, 
  public sphlatch::SPHfluidGhost {};

typedef particle                               partT;
typedef ghost                                  ghstT;


#include "bhtree_treedump.cpp"
typedef sphlatch::BHTreeDump                   dumpT;

#include "bhtree_worker_grav.cpp"
typedef sphlatch::fixThetaMAC                  macT;
typedef sphlatch::GravityWorker<macT, partT>   gravT;

//#include "bhtree_worker_sphsum.cpp"

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

   const size_t       noParts = 1000000;
   std::vector<partT> particles(noParts);
   for (size_t i = 0; i < noParts; i++)
   {
      particles[i].pos[0] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[1] = static_cast<fType>(rand()) / RAND_MAX;
      particles[i].pos[2] = static_cast<fType>(rand()) / RAND_MAX;

      particles[i].m = 1.;

      particles[i].id   = i;
      particles[i].cost = 1.;

      Tree.insertPart( particles[i] );
   }

   //Tree.insertParts(particles);
   std::cout << "Tree.update() .........\n";
   Tree.update();
   std::cout << "Tree.update() finished.\n";

   //dumper.ptrDump("postCZdump.txt");
   //dumper.dotDump("postCZdump.dot");

   /*testWorker.dispRoot();
      testWorker.build();
      testWorker.test();*/

   //dumper.dotDump("pre__move.dot");
   //dumper.ptrDump("pre__move.txt");
   for (size_t i = 0; i < noParts; i++)
   {
      if (particles[i].pos[2] < 0.5)
         particles[i].pos[2] += 0.5;

      particles[i].cost = 1.;
   }

   std::cout << "Tree.update() .........\n";
   Tree.update();
   std::cout << "Tree.update() finished.\n";

   //std::cout << (static_cast<std::vector<sphlatch::treeGhost> >(particles))[0].pos << "\n";

   //dumper.dotDump("post_move.dot");
   //dumper.ptrDump("post_move.txt");

   gravT gravWorker(&Tree);

   std::cout << "gravity ...\n";
   treeT::czllPtrVectT CZbottomLoc = Tree.getCZbottomLoc();
   const int noCZbottomLoc = CZbottomLoc.size();

   omp_set_num_threads(2);
//#pragma omp parallel for firstprivate(gravWorker)
   for (int i = 0; i < noCZbottomLoc; i++)
   {
     //std::cout << "  grav: " << CZbottomLoc[i] << "\n";
     gravWorker.calcGravity(CZbottomLoc[i]);
   }
   std::cout << "gravity finished \n";
   
   //for (size_t i = 0; i < noParts; i++)
     //std::cout << particles[i].acc << "\n";

   MPI::Finalize();
   return(0);
}

#include <iostream>
#include <vector>

#include <omp.h>

#include "typedefs.h"
typedef sphlatch::fType                  fType;
typedef sphlatch::vect3dT                vect3dT;

#include "bhtree_particle.h"
typedef sphlatch::treeGhost              treeGhostT;
typedef sphlatch::treePart               treePartT;

#include "sph_fluid_particle.h"
typedef sphlatch::SPHfluidGhost          sphGhostT;
typedef sphlatch::SPHfluidPart           sphPartT;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager   commT;

#include "bhtree_node_particle.h"

class ghostT : public treeGhostT, public sphGhostT
{
public:
};

class partT : public treePartT, public sphPartT
{
public:
  void swap(partT& _swpPart)
  {
    static partT _tmp; 
    _tmp = _swpPart;
    _swpPart = *this;
    *this = _tmp;

    _swpPart.treeNode->partPtr = &_swpPart;
    treeNode->partPtr = this;
  }
};

#include "particle_allocator.h"


int main(int argc, char* argv[])
{
   MPI::Init(argc, argv);

   //treeT& Tree(treeT::instance());
   commT& Comm(commT::instance());

   ghostT ghost;
   partT  part;

   partT* partPtr = &part;
   
   treePartT* treePartPtr = &part;

   commT::DataTypeCreator creator(&ghost);
   creator += &ghost.pos;
   creator += &ghost.m;
   commT::dataType ghostPos = creator.finalize();

   creator(&part);
   creator += &part.m;
   creator += &part.pos;
   commT::dataType partPos = creator.finalize();

   const int rank = MPI::COMM_WORLD.Get_rank();

   if (rank == 0)
   {
      part.m   = 3.14e98;
      part.pos = 1, 2, 3;

      MPI::COMM_WORLD.Send(&part, 1, partPos, 1, 0);
   }
   else if (rank == 1)
   {
      ghost.m   = 3.;
      ghost.pos = 0, 0, 0;
      
      std::cout << "rank" << rank
                << ": " << ghost.m
                << "   pos" << ghost.pos
                << "\n";

      MPI::COMM_WORLD.Recv(&ghost, 1, ghostPos, 0, 0);
      std::cout << "rank" << rank
                << ": " << ghost.m
                << "   pos" << ghost.pos
                << "\n";
   }

   MPI::Finalize();
   return(0);
}

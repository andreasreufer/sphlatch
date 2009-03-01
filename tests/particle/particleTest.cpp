#include <iostream>
#include <vector>

#include <omp.h>

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

#include "bhtree_particle.h"
typedef sphlatch::treeGhost           partT;

#include "communication_manager.h"
typedef sphlatch::CommunicationManager   commT;

int main(int argc, char* argv[])
{
   MPI::Init(argc, argv);

   //treeT& Tree(treeT::instance());
   commT& Comm(commT::instance());

   partT part;

   //sphlatch::SPHpart sph_part;

   const MPI::Aint addr_base = MPI::Get_address(&part);
   const MPI::Aint addr_mass = MPI::Get_address(&part.m);

   const MPI::Aint offs_mass = addr_mass - addr_base;

   const int blength[3] = {1,1,1};
   const MPI::Aint displac[3] = {0, offs_mass, offs_mass+sizeof(fType)};

   const MPI::Datatype dtypes[3] = {MPI::LB, MPI::DOUBLE, MPI::UB};
   
   static MPI::Datatype MPImassT = MPI::Datatype::Create_struct(3, blength, displac, dtypes);

   MPImassT.Commit();

   //MPI::Datatype::Commit(MPImassT);

   const int rank = MPI::COMM_WORLD.Get_rank();
   
   if ( rank == 0 )
   {
     part.m = 3.14e98;
     part.pos = 1,2,3;

     MPI::COMM_WORLD.Send( &part, 1, MPImassT, 1, 0);
   }
   else if ( rank == 1 )
   {
     part.m = 0.;
     part.pos = 0,0,0;
     
     MPI::COMM_WORLD.Recv( &part, 1, MPImassT, 0, 0);
   }

   std::cout << "rank" << rank << ": " << part.m << "   pos" << part.pos << "\n";


   MPI::Finalize();
   return(0);
}

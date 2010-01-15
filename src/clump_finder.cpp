#ifndef SPHLATCH_CLUMP_FINDER_CPP
#define SPHLATCH_CLUMP_FINDER_CPP

/*
 *  clump_finder.cpp
 *
 *
 *  Created by Andreas Reufer on 10.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "particle_set.cpp"
#include "clump_particle.h"

#include "bhtree.cpp"
typedef sphlatch::BHTree   treeT;

#include "bhtree_worker_neighfunc.cpp"

namespace sphlatch {

template<typename _partT>
void getClumps(ParticleSet<_partT>& _parts)
{
   treeT& Tree(treeT::instance());

   const size_t nop       = _parts.getNop();
   const fType  costppart = 1. / nop;

   Tree.setExtent(_parts.getBox() * 1.1);

   for (size_t i = 0; i < nop; i++)
   {
      _parts[i].cost = costppart;
      Tree.insertPart(_parts[i]);
   }
   Tree.update(0.8, 1.2);

   cType nextFreeId = nop + 1;

   NeighFindWorker<_partT> NFW(&Tree);
   for (size_t i = 0; i < nop; i++)
   {
       // neigh search & add to friends stack

     /*switch(_parts[i].clumpid)
     {
       case CLUMPNONE:
         break;

       case CLUMPNOTSET:
       // get new clump id
       // neigh search & add to friends stack
       // infect friends
         break;
       
       default:
       // rho > rhoMin?
       // neigh search & add to friends stack
       // infect friends
         break;

     }*/

     const cType clumpid = _parts[i].clumpid;

     if ( clumpid != CLUMPNONE )
     {
       cType curid;
       if ( clumpid != CLUMPNOTSET )
       {
         curid = nextFreeId;
         nextFreeId++;
       }
       else
         curid = clumpid;

       //std::list<_partT*> pFriends = NF(&(_parts[i]), _parts[i].h);

       std::list<_partT*> ffriends;
       ffriends.push_back(&(_parts[i]));

       while ( ffriends.size() != 0 )
       {
         const _partT* curff = ffriends.back();
         ffriends.pop_back();

         // condition for a friend
         //if curff->rho > rhoMin;
       }


     }
   }

   //treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   //const int           noCZbottomLoc = CZbottomLoc.size();

   Tree.clear();
}


/*class clump
   {
   public:
   const vect3dT pos, vel, acc, L;
   };*/


/*   void flavourClumps()
   { }

   void findClumps()
   { }
   }*/
}

#endif

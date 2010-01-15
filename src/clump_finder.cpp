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
void getClumps(ParticleSet<_partT>& _parts, const fType _rhoMin)
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
      std::cout << "particle " << i << " rho = " << _parts[i].rho << "\n";
      cType& clumpid(_parts[i].clumpid);

      //
      // if density is too low, particle does not belong to any clump
      // regardless of previous state
      //
      if (_parts[i].rho < _rhoMin)
      {
         clumpid = CLUMPNONE;
         continue;
      }

      //
      // density is high enough to belong to a clump
      //
      if ((clumpid == CLUMPNONE) || (clumpid == CLUMPNOTSET))
      {
         clumpid = nextFreeId;
         nextFreeId++;
      }

      //
      // now start to look for friends
      //
      std::list<_partT*> ffriends;
      ffriends.push_back(&(_parts[i]));

      const cType curclump = clumpid;
      while (ffriends.size() != 0)
      {
         std::cout << " ff: " << ffriends.size() << "\n";
         _partT* const curff = ffriends.back();
         ffriends.pop_back();
         std::cout << " ff: " << ffriends.size() << "\n";


         if (curff->rho > _rhoMin)
         {
            if (curff->clumpid > 0)
            {
               // add to dictionnary
            }
            curff->clumpid = curclump;
            // add neighbours to the friends of friends list
            std::list<_partT*> newffriends = NFW(curff, curff->h);
            std::cout << "   " << newffriends.size() << " new friends found\n";

            // only add the particle, if it does not yet belong to a clump
            //ffriends.splice(ffriends.end(), newffriends);
         }

         // condition for a friend
         //if curff->rho > rhoMin;
      }
   }

   // assign real clump IDs

   // output clumps?

   // clear tree
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

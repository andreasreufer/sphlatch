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

#include "bhtree.cpp"
typedef sphlatch::BHTree   treeT;

#include "bhtree_worker_neighfunc.cpp"

namespace sphlatch {
enum clumpType
{
   CLUMPNONE = -1, CLUMPNOTSET = 0
};

class clumpPart
{
   iType clumpID;
};

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


   //treeT::czllPtrVectT CZbottomLoc   = Tree.getCZbottomLoc();
   //const int           noCZbottomLoc = CZbottomLoc.size();

   Tree.clear();
}


/*class clump
   {
   public:
   const vect3dT pos, vel, acc, L;
   };*/

/*template<typename _partT>
   struct clumpFlavourer
   {

   void   preSum(_partT* const _i)
   {
      rhoi = 0.;
   }

   void operator()(_partT* const _i,
                   const _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      const fType r   = sqrt(_rr);
   }

   void postSum(_partT* const _i)
   {
      _i->rho = rhoi;
   }

   void flavourClumps()
   { }

   void findClumps()
   { }
   }*/
}

#endif

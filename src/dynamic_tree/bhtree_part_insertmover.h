#ifndef BHTREE_PART_INSERTMOVER_H
#define BHTREE_PART_INSERTMOVER_H

/*
 *  bhtree_part_insertmover.h
 *
 *  Created by Andreas Reufer on 02.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"
#include "bhtree_errhandler_tmp.h"

namespace sphlatch {
class BHTreePartsInsertMover : public BHTreeWorker {
public:

   typedef BHTreePartsInsertMover   selfT;

   BHTreePartsInsertMover(treePtrT _treePtr) 
     : BHTreeWorker(_treePtr) { }
   BHTreePartsInsertMover(const selfT& _inserter) 
     : BHTreeWorker(_inserter) { }
   ~BHTreePartsInsertMover() { }

public:
   void insert(const size_t _i);
   void moveAll();

private:
   void recursor(const size_t _i);
};

///
/// entry function to insert particle:
/// - check whether particle lies inside root cell
/// - if so, call insertion recursor
///
void BHTreePartsInsertMover::insert(const size_t _i)
{
   goRoot();

   ///
   /// check whether the particles comes to lie in the root
   /// cell, if not throw an exception
   ///
   const fType curX = pos(_i, X);
   const fType curY = pos(_i, Y);
   const fType curZ = pos(_i, Z);
   const fType curM = m(_i);

   if (not pointInsideCell(curX, curY, curZ))
      throw PartOutsideTree(_i, rootPtr);

   ///
   /// allocate particle node, wire it to its proxy and
   /// set nodes parameters
   ///
   const partPtrT curPartPtr = treePtr->partAllocator.pop();
   treePtr->partProxies[_i] = curPartPtr;

   static_cast<partPtrT>(curPartPtr)->clear();
   static_cast<partPtrT>(curPartPtr)->xPos = curX;
   static_cast<partPtrT>(curPartPtr)->yPos = curY;
   static_cast<partPtrT>(curPartPtr)->zPos = curZ;
   static_cast<partPtrT>(curPartPtr)->mass = curM;

   static_cast<partPtrT>(curPartPtr)->ident = _i;
   static_cast<partPtrT>(curPartPtr)->depth = 0;

   ///
   /// wire the root cell as particles parent (in case the tree
   /// is not cleaned before use, this particle is not found in
   /// a tree walk, as it is not yet wired as a child of a cell)
   ///
   static_cast<partPtrT>(curPartPtr)->parent = treePtr->rootPtr;
}

///
/// particle insertion recursor
///
///

/*//void BHTreePartsInsertMover::recursor(const size_t _i)
   void BHTreePartsInsertMover::recursor(const partPtrT _ptr)
   {
   const fType curX    = pos(_i, X);
   const fType curY    = pos(_i, Y);
   const fType curZ    = pos(_i, Z);
   const fType curM    = m(_i);
   const fType curCost = cost(_i);

   const size_t tarOct = getOctant(curX, curY, curZ);

   static_cast<cellPtrT>(curPtr)->cost += curCost;
   //static_cast<cellPtrT>(curPtr)->nPar++;

   ///
   /// If targeted child is empty, place the particle there
   ///
   if (static_cast<cellPtrT>(curPtr)->child[tarOct] == NULL)
   {
      ///
      /// wire new particle node, go down to it and clear
      /// it (you never know who touched it before)
      ///
      const partPtrT partPvt = treePtr->partAllocator.pop();
      partPvt->parent = curPtr;
      static_cast<cellPtrT>(curPtr)->child[tarOct] = partPvt;

      goChild(tarOct);

      static_cast<partPtrT>(curPtr)->clear();

      treePtr->partProxies[_i] = static_cast<partPtrT>(curPtr);

      static_cast<partPtrT>(curPtr)->ident = _i;
      static_cast<partPtrT>(curPtr)->depth = curPtr->parent->depth + 1;


      goUp();
   }

   ///
   /// ... or if existing child is a node, then try to place the particle
   /// as a child of this node
   ///
   else if (not static_cast<cellPtrT>(curPtr)->child[tarOct]->isParticle)
   {
      goChild(tarOct);
      recursor(_i);
      goUp();
   }

   ///
   /// ... or if the existing child is a particle, then replace it by a
   /// cell node and try to insert them both again at the current node
   ///
   else
   {
      ///
      /// goto child, delete old particle but save its index
      /// for later insertion
      ///
      goChild(tarOct);
      const identType old_i = static_cast<partPtrT>(curPtr)->ident;
      goUp();
      treePtr->partAllocator.push(
         static_cast<partPtrT>(
            static_cast<cellPtrT>(curPtr)->child[tarOct]));

      ///
      /// new cell
      ///

      ///
      /// try inserting particles again
      ///
      goChild(tarOct);

      recursor(old_i);
      recursor(_i);
      goUp();

 */

/*
   goChild(targetOctant);
   size_t residentIdx = static_cast<partPtrT>(curNodePtr)->ident;
   goUp();

   /// replace particle by cell node
   delete static_cast<cellPtrT>(curNodePtr)->child[targetOctant];
   static_cast<cellPtrT>(curNodePtr)->child[targetOctant] = NULL;
   newCellChild(targetOctant);

   /// and try to insert both particles again
   goChild(targetOctant);
   ///
   /// check whether we are too deep
   ///
   if (curNodePtr->depth > 128)
   throw PartsTooClose(curNodePtr->depth,
                      residentIdx,
                      _newPartIdx,
                      rootPtr);

   insertParticleRecursor(residentIdx);
   insertParticleRecursor(_newPartIdx);
   goUp();*/
//}
//}
};

#endif

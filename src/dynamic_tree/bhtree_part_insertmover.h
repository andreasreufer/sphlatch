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
   void pushToCZ();
   void pushDown();

private:
   void pushToCZsingle(const partPtrT _partPtr);
   void pushDownSingle(const partPtrT _partPtr);

   //void recursor(const size_t _i);
   //void move(const size_t _i);
};

///
/// entry function to insert particle:
/// - check whether particle lies inside root cell
/// - if so, call insertion recursor
///
void BHTreePartsInsertMover::insert(const size_t _i)
{
   ///
   /// allocate particle node, wire it to its proxy and
   /// set nodes parameters
   ///
   const partPtrT newPart = treePtr->partAllocator.pop();

   treePtr->partProxies[_i] = newPart;

   newPart->ident     = _i;
   newPart->depth     = 0;
   newPart->isSettled = false;

   newPart->xPos = pos(_i, X);
   newPart->yPos = pos(_i, Y);
   newPart->zPos = pos(_i, Z);
   newPart->mass = m(_i);

   ///
   /// wire the root cell as particles parent (in case the tree
   /// is not cleaned before use, this particle is not found in
   /// a tree walk, as it is not yet wired as a child of a cell)
   ///
   newPart->parent = treePtr->rootPtr;
}

///
/// update particle position and move it to the CZ bottom
///
void BHTreePartsInsertMover::pushToCZsingle(const partPtrT _partPtr)
{
   ///
   /// get particle index, its child number and check whether
   /// with the new position the particle still lies in the same
   /// cell. if yes, do nothing and return.
   ///
   curPtr = _partPtr->parent;
   const size_t oldOct = getOctant(_partPtr);

   assert(oldOct >= 0 && oldOct < 8);

   const size_t partIdx = _partPtr->ident;
   const fType  posX    = pos(partIdx, X);
   const fType  posY    = pos(partIdx, Y);
   const fType  posZ    = pos(partIdx, Z);
   const fType  mass    = m(partIdx);

   _partPtr->xPos = posX;
   _partPtr->yPos = posY;
   _partPtr->zPos = posZ;
   _partPtr->mass = mass;

   if (pointInsideCell(posX, posY, posZ) &&
       (getOctant(posX, posY, posZ) == oldOct))
      return;

   ///
   /// particles does not lie in the same cell as beforce, so
   /// start to move it
   ///
   static_cast<gcllPtrT>(curPtr)->child[oldOct] = NULL;

   ///
   /// check whether the particle lies inside the root cell
   ///
   if (not pointInsideCell(posX, posY, posZ, rootPtr))
   { }

   ///
   /// go up until the particle lies in the current cell. as
   /// soon as we hit a CZ cell, start substracting relative
   /// cost from each CZ cell we encounter
   ///
   bool  CZencounter = false;
   fType partCost    = 0.;
   while (not pointInsideCell(posX, posY, posZ))
   {
      if (not CZencounter)
      {
         if (curPtr->isCZ)
         {
            CZencounter = true;
            partCost    =
               static_cast<czllPtrT>(curPtr)->cost /
               static_cast<czllPtrT>(curPtr)->noParts;
         }
      }
      else
      {
         /// use of relativ cost necessary
         static_cast<czllPtrT>(curPtr)->noParts--;
         static_cast<czllPtrT>(curPtr)->cost -= partCost;
      }

      goUp();
   }

   ///
   /// in case are now in the CZ tree, try to go down again to
   /// the bottom layer
   ///
   if (CZencounter)
   {
      while (not curPtr->atBottom)
         goChild(getOctant(posX, posY, posZ));
   }

   _partPtr->parent = curPtr;
}

///
/// update particle position and move it to the CZ bottom
///
void BHTreePartsInsertMover::pushDownSingle(const partPtrT _partPtr)
{
   curPtr = _partPtr->parent;

   const fType posX = _partPtr->xPos;
   const fType posY = _partPtr->yPos;
   const fType posZ = _partPtr->zPos;

   assert(pointInsideCell(posX, posY, posZ));
   assert(not curPtr->isParticle);

   bool isSettled = false;
   while (not isSettled)
   {
      const size_t curOct = getOctant(posX, posY, posZ);
      ///
      /// child is empty, particle can be inserted directly
      ///
      if (static_cast<cellPtrT>(curPtr)->child[curOct] == NULL)
      {
         static_cast<cellPtrT>(curPtr)->child[curOct] = _partPtr;

         _partPtr->parent    = curPtr;
         _partPtr->depth     = curPtr->depth + 1;
         _partPtr->isSettled = true;

         isSettled = true;
      }
      else
      {
         const nodePtrT curParentPtr = curPtr;
         goChild(curOct);

         ///
         /// 
         ///
         if (static_cast<cellPtrT>(curPtr)->child[curOct]->isParticle)
         {
            const partPtrT resPartPtr = static_cast<partPtrT>(curPtr);

            curPtr = treePtr->cellAllocator.pop();
            static_cast<cellPtrT>(curPtr)->clear();
            
            curPtr->parent = curParentPtr;
            static_cast<cellPtrT>(curPtr)->inheritCellPos(curOct);

            pushDownSingle(resPartPtr);
         }
      }
   }
}

};

#endif

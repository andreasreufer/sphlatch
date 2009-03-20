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
//#include "bhtree_errhandler_tmp.h"

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
   void insert(partT& _part);

   //void insert(std::vector<particle>& _parts);
   //void insert(std::vector<ghost>& _ghosts);
   void pushToCZ();

   // temporary
   void pushDown(const partT& _part);

private:
   void pushToCZsingle(const pnodPtrT _pnodPtr);
   void pushUpAndToCZsingle(const pnodPtrT _pnodPtr);
   void pushDownSingle(const pnodPtrT _pnodPtr);
};

///
/// entry function to insert particle:
/// - check whether particle lies inside root cell
/// - if so, call insertion recursor
///
void BHTreePartsInsertMover::insert(partT& _part)
{
   ///
   /// allocate particle node, wire it to its proxy and
   /// set nodes parameters
   ///
   const pnodPtrT newPartPtr = new pnodT;

   newPartPtr->clear();

   newPartPtr->partPtr = &_part;
   _part.treeNode      = newPartPtr;
   newPartPtr->update();

   newPartPtr->ident     = _part.id;
   newPartPtr->depth     = 0;
   newPartPtr->isSettled = false;

   ///
   /// wire the root cell as particles parent (in case the tree
   /// is not cleaned before use, this particle is not found in
   /// a tree walk, as it is not yet wired as a child of a cell)
   ///
   newPartPtr->parent = treePtr->rootPtr;
   static_cast<gcllPtrT>(newPartPtr->parent)->noParts++;

   //pushUpAndToCZsingle(newPartPtr);
}

///
/// update particle position and move it to the current CZ bottom
///
void BHTreePartsInsertMover::pushUpAndToCZsingle(const pnodPtrT _pnodPtr)
{
   ///
   /// get particle index, its child number and check whether
   /// with the new position the particle still lies in the same
   /// cell. if yes, do nothing and return.
   ///
   curPtr = _pnodPtr->parent;
   const size_t oldOct = getChildNo(_pnodPtr);

   assert(oldOct >= 0 && oldOct < 8);

   ///
   /// update particle position and mass
   ///
   _pnodPtr->update();
   const vect3dT pos = _pnodPtr->pos;

   if (pointInsideCell(pos) &&
       (getOctant(pos) == oldOct))
      return;

   ///
   /// particles does not lie in the same cell as beforce, so
   /// start to move it
   ///
   static_cast<gcllPtrT>(curPtr)->child[oldOct] = NULL;

   ///
   /// check whether the particle lies inside the root cell
   ///
   if (not pointInsideCell(pos, rootPtr))
   { }

   ///
   /// go up until the particle lies in the current cell. as
   /// soon as we hit a CZ cell, start substracting cost from
   /// each CZ cell we encounter
   ///
   bool        CZencounter = false;
   const fType partCost    = static_cast<pnodPtrT>(_pnodPtr)->partPtr->cost;
   while (not pointInsideCell(pos))
   {
      if (not CZencounter)
      {
         if (curPtr->isCZ)
            CZencounter = true;
      }

      static_cast<gcllPtrT>(curPtr)->cost -= partCost;
      static_cast<gcllPtrT>(curPtr)->noParts--;

      goUp();
   }

   ///
   /// in case are now in the CZ tree, try to go down again to
   /// the bottom layer
   ///
   if (CZencounter)
   {
      while (not curPtr->atBottom)
      {
         goChild(getOctant(pos));
         static_cast<gcllPtrT>(curPtr)->cost += partCost;
         static_cast<gcllPtrT>(curPtr)->noParts++;
      }
   }

   _pnodPtr->parent = curPtr;
}

///
/// push down a single particle
///
void BHTreePartsInsertMover::pushDown(const partT& _part)
{
   pushDownSingle(_part.treeNode);
}

void BHTreePartsInsertMover::pushDownSingle(const pnodPtrT _pnodPtr)
{
   curPtr = _pnodPtr->parent;
   const vect3dT pos = _pnodPtr->pos;
   const fType partCost = _pnodPtr->partPtr->cost;

   assert(pointInsideCell(pos));
   assert(not curPtr->isParticle);

   bool isSettled = false;
   while (not isSettled)
   {
      const size_t curOct = getOctant(pos);
      ///
      /// child is empty, particle can be inserted directly
      ///
      if (static_cast<cellPtrT>(curPtr)->child[curOct] == NULL)
      {
         static_cast<cellPtrT>(curPtr)->child[curOct] = _pnodPtr;

         _pnodPtr->parent    = curPtr;
         _pnodPtr->depth     = curPtr->depth + 1;
         _pnodPtr->isSettled = true;

         isSettled = true;
      }
      else
      {
         if (static_cast<cellPtrT>(curPtr)->child[curOct]->isParticle)
            partToCell(curPtr, curOct);

         goChild(curOct);
         
         static_cast<gcllPtrT>(curPtr)->cost += partCost;
         static_cast<gcllPtrT>(curPtr)->noParts++;
      }
   }
}
};

#endif

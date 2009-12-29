#ifndef BHTREE_PART_INSERTMOVER_CPP
#define BHTREE_PART_INSERTMOVER_CPP

/*
 *  bhtree_part_insertmover.cpp
 *
 *  Created by Andreas Reufer on 02.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_part_insertmover.h"
#include "bhtree_worker.cpp"

namespace sphlatch {
BHTreePartsInsertMover::BHTreePartsInsertMover(treePtrT _treePtr) :
   BHTreeWorker(_treePtr),
   maxDepth(_treePtr->maxDepth)
{ }

BHTreePartsInsertMover::BHTreePartsInsertMover(const selfT& _inserter) :
   BHTreeWorker(_inserter),
   maxDepth(_inserter.maxDepth)
{ }

BHTreePartsInsertMover::~BHTreePartsInsertMover() { }


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

   newPartPtr->ident     = _part.id;
   newPartPtr->depth     = 0;
   newPartPtr->isSettled = false;

   ///
   /// wire the root cell as particles parent (in case the tree
   /// is not cleaned before use, this particle is not found in
   /// a tree walk, as it is not yet wired as a child of a cell)
   ///
   newPartPtr->parent = treePtr->rootPtr;

   static_cast<czllPtrT>(rootPtr)->noParts++;
   static_cast<czllPtrT>(rootPtr)->relCost += _part.cost;

   pushUpAndToCZSingle(newPartPtr);
}

void BHTreePartsInsertMover::move(const czllPtrT _czll)
{
   nodePtrT curPart = _czll->chldFrst;

   if (_czll->chldLast != NULL)
      _czll->chldLast->next = NULL;

   ///
   /// if there is a last child, set its next pointer to NULL,
   /// so that the next walk for moving terminates
   ///
   while (curPart != NULL)
   {
      const nodePtrT nextPart = curPart->next;
      if (curPart->isParticle)
         pushUpAndToCZSingle(static_cast<pnodPtrT>(curPart));
      curPart = nextPart;
   }
}

void BHTreePartsInsertMover::pushDownOrphans(const czllPtrT _czll)
{
   nodePtrT curPart = _czll->orphFrst;

   while (curPart != NULL)
   {
      pushDownSingle(static_cast<pnodPtrT>(curPart));
      curPart = curPart->next;
   }

   _czll->orphFrst = NULL;
   _czll->orphLast = NULL;
}

///
/// update particle position and move it to the current CZ bottom
///
void BHTreePartsInsertMover::pushUpAndToCZSingle(const pnodPtrT _pnodPtr)
{
   ///
   /// get particle index, its child number and check whether
   /// with the new position the particle still lies in the same
   /// cell. if yes, do nothing and return.
   ///
   curPtr = _pnodPtr->parent;
   const size_t oldOct = getChildNo(_pnodPtr);

   ///
   /// update particle position and mass
   ///
   _pnodPtr->update();
   const vect3dT pos = _pnodPtr->pos;

   ///
   /// a short-cut for those particles which
   /// will stay in the same cell octant
   ///
   if (pointInsideCell(pos) &&
       (getOctant(pos) == oldOct))
      return;

   ///
   /// if the particle was a child of the parent cell
   ///
   if (oldOct != 8)
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
      CZencounter |= curPtr->isCZ;

      if (CZencounter)
      {
         static_cast<czllPtrT>(curPtr)->relCost -= partCost;
         static_cast<czllPtrT>(curPtr)->noParts--;
      }
      goUp();
   }
   _pnodPtr->parent = curPtr;

   ///
   /// in case are now in the CZ tree, try to go down again to
   /// the bottom layer
   ///
   if (CZencounter)
   {
      while (not curPtr->atBottom)
      {
         goChild(getOctant(pos));
         static_cast<czllPtrT>(curPtr)->relCost += partCost;
         static_cast<czllPtrT>(curPtr)->noParts++;
      }
      static_cast<czllPtrT>(curPtr)->adopt(_pnodPtr);
   }
   else
      pushDownSingle(_pnodPtr);
}

void BHTreePartsInsertMover::pushDownSingle(const pnodPtrT _pnodPtr)
{
   curPtr = _pnodPtr->parent;
   const vect3dT pos = _pnodPtr->pos;

   //if (not pointInsideCell(pos))
   //  std::cout << _pnodPtr->partPtr->id << "\n";

   assert(pointInsideCell(pos));
   assert(not curPtr->isParticle);

   bool isSettled = false;
   while (not isSettled)
   {
      const size_t curOct = getOctant(pos);
      ///
      /// child is empty, particle can be inserted directly
      ///
      if (static_cast<gcllPtrT>(curPtr)->child[curOct] == NULL)
      {
         static_cast<gcllPtrT>(curPtr)->child[curOct] = _pnodPtr;

         _pnodPtr->parent    = curPtr;
         _pnodPtr->depth     = curPtr->depth + 1;
         _pnodPtr->isSettled = true;

         isSettled = true;
      }
      else
      {
         /*if (curPtr->depth > maxDepth)
         {
            std::cerr
             << static_cast<pnodPtrT>(static_cast<gcllPtrT>(curPtr)->child[
                                         curOct])->pos << "\n"
                                                       << pos << "\n"
                                                       << static_cast<pnodPtrT>(
               static_cast<gcllPtrT>(curPtr)->child[
                  curOct
               ])->pos - pos << "\n"
                             <<
            static_cast<gcllPtrT>(curPtr)->cen << " " <<
            static_cast<gcllPtrT>(curPtr)->clSz << "\n";
            exit(1);
         }*/
         if (static_cast<gcllPtrT>(curPtr)->child[curOct]->isParticle)
            static_cast<gcllPtrT>(curPtr)->child[curOct] =
               partToCell(curPtr, curOct);
         goChild(curOct);
      }
   }
}
};

#endif

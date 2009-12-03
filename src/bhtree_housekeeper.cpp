#ifndef BHTREE_HOUSEKEEPER_CPP
#define BHTREE_HOUSEKEEPER_CPP

/*
 *  bhtree_housekeeper.cpp
 *
 *  Created by Andreas Reufer on 14.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <vector>

#include "bhtree_housekeeper.h"
#include "bhtree_worker.cpp"

namespace sphlatch {
BHTreeHousekeeper::BHTreeHousekeeper(const treePtrT _treePtr) :
   BHTreeWorker(_treePtr)
{ }

BHTreeHousekeeper::BHTreeHousekeeper(const BHTreeHousekeeper& _hk) :
   BHTreeWorker(_hk)
{ }

BHTreeHousekeeper::~BHTreeHousekeeper()
{ }

void BHTreeHousekeeper::setNext(const czllPtrT _czll)
{
   ///
   /// wire next pointer by doing a preorder tree walk
   ///
   curPtr  = _czll;
   lastPtr = _czll;
   setNextRecursor();
   lastPtr->next = NULL;

   ///
   /// wire the child first and child last pointers
   ///
   _czll->chldFrst = _czll->next;
   _czll->chldLast = lastPtr;

   curPtr = _czll;
}

void BHTreeHousekeeper::setNextCZ()
{
   if ( rootPtr->atBottom )
     return;

   curPtr  = rootPtr;
   lastPtr = rootPtr;
   setNextCZRecursor();
   lastPtr->next = NULL;
}

void BHTreeHousekeeper::setSkip()
{
   ///
   /// do a preorder walk by using the next pointers
   /// and store at each depth the last cell encountered
   /// in a list
   ///
   /// when going up again to depth n and encountering a cell,
   /// let the "skip"-pointer of each cell in the list with
   /// >= n point at the current cell
   ///

   ///
   /// prepare the lastSkipee list
   ///
   const size_t maxDepth = treePtr->maxDepth;

   if (lastSkipeeAtDepth.size() != maxDepth)
      lastSkipeeAtDepth.resize(maxDepth);
   for (size_t i = 0; i < maxDepth; i++)
      lastSkipeeAtDepth[i] = NULL;

   ///
   /// do the next-walk
   ///
   goRoot();
   goNext();
   while (curPtr != NULL)
   {
      const size_t depth = curPtr->depth;

      size_t i = depth;
      while (lastSkipeeAtDepth[i] != NULL)
      {
         lastSkipeeAtDepth[i]->skip = static_cast<gcllPtrT>(curPtr);
         lastSkipeeAtDepth[i]       = NULL;
         i++;
      }
      if (not curPtr->isParticle)
         lastSkipeeAtDepth[depth] = static_cast<gcllPtrT>(curPtr);

      goNext();
   }
}

void BHTreeHousekeeper::setNextRecursor()
{
   lastPtr->next = curPtr;
   lastPtr       = curPtr;

   if (not curPtr->isParticle)
   {
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            goChild(i);
            setNextRecursor();
            goUp();
         }
      }
   }
}

void BHTreeHousekeeper::setNextCZRecursor()
{
   if (not curPtr->isParticle)
   {
      lastPtr->next = curPtr;

      if (curPtr->atBottom)
      {
         lastPtr = static_cast<czllPtrT>(curPtr)->chldLast;
      }
      else
      {
         lastPtr = curPtr;

         for (size_t i = 0; i < 8; i++)
         {
            if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
            {
               goChild(i);
               setNextCZRecursor();
               goUp();
            }
         }
      }
   }
}

void BHTreeHousekeeper::minTree(const czllPtrT _czll)
{
   if (_czll->chldFrst == NULL)
      return;

   curPtr = _czll;

   const nodePtrT chldLastNext = _czll->chldLast->next;

   nodePtrT lastOk = NULL, nextOk = NULL, nextChld;

   // maybe the last node is deleted during housekeeping, so
   // we introduce a dummy cell won't be deleted
   const nodePtrT chldLastDummy = new pnodT;
   _czll->chldLast->next = chldLastDummy;

   while (curPtr != chldLastDummy)
   {
      nextChld = curPtr->next;

      if (nextChld->isParticle || nextChld->isCZ)
      {
         lastOk = curPtr;
         goNext();
         continue;
      }

      const size_t noChld = static_cast<gcllPtrT>(nextChld)->getNoChld();

      switch (noChld)
      {
      case 0:
      {
         // the cell has no child, so delete it directly
         nextOk = nextChld->next;
         const size_t wasChild = getChildNo(nextChld, nextChld->parent);
         static_cast<gcllPtrT>(nextChld->parent)->child[wasChild] = NULL;
         delete static_cast<qcllPtrT>(nextChld);

         curPtr->next = nextOk;
         break;
      }

      case 1:
      {
         // only delete chains not leading anywhere
         nodePtrT chainee = nextChld;

         while (not chainee->isParticle &&
                not chainee->isCZ &&
                static_cast<gcllPtrT>(chainee)->getNoChld() == 1)
            chainee = chainee->next;

         if (chainee->isParticle || chainee->isCZ ||
             (not static_cast<gcllPtrT>(chainee)->getNoChld() == 0))
         {
            curPtr = chainee;
            break;
         }

         const size_t wasChild = getChildNo(nextChld, nextChld->parent);
         static_cast<gcllPtrT>(nextChld->parent)->child[wasChild] = NULL;

         // delete the chainees, until we reach the last chainee
         while (nextChld != chainee)
         {
            nextOk = nextChld->next;
            delete nextChld;
            nextChld = nextOk;
         }

         // set next ptr of the current node to the next ptr of the last
         // chainee
         curPtr->next = nextChld->next;
         delete chainee;

         break;
      }

      default:
      {
         goNext();
         break;
      }
      }
   }


   delete chldLastDummy;
   lastOk->next    = chldLastNext;
   _czll->chldLast = lastOk;
}
};

#endif

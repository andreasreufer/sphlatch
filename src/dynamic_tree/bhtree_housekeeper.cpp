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
      if (not curPtr->isParticle)
      {
         size_t i = depth;
         while (lastSkipeeAtDepth[i] != NULL)
         {
            lastSkipeeAtDepth[i]->skip = static_cast<gcllPtrT>(curPtr);
            lastSkipeeAtDepth[i]       = NULL;
            i++;
         }
         lastSkipeeAtDepth[depth] = static_cast<gcllPtrT>(curPtr);
      }
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
   std::cout << "\n\n\nminimize czll " << _czll << "\n";
   
   if (_czll->chldFrst == NULL)
      return;
   //curPtr = _czll->chldFrst;
   curPtr = _czll;
   
   const nodePtrT chldLastNext = _czll->chldLast->next;
   const nodePtrT chldLast = _czll->chldLast;

   std::cout << "chldLast: " << chldLast << " " << chldLast->next << "\n";

   nodePtrT lastOk = NULL, nextOk = NULL, nextChld;

   // maybe the last node is deleted during housekeeping, so
   // we introduce a dummy cell won't be deleted
   const nodePtrT chldLastDummy = new pnodT;
   _czll->chldLast->next = chldLastDummy;


   while (curPtr != chldLastDummy)
   {
      nextChld = curPtr->next;
      std::cout << __LINE__ << " beg  " << curPtr << " " << nextChld << "\n";

      if (nextChld->isParticle || nextChld->isCZ)
      {
         std::cout << __LINE__ << " skip " << nextChld << " " << nextChld->isParticle << " " << nextChld->isCZ << "\n";
         lastOk = curPtr;
         goNext();
         continue;
      }

      const size_t noChld = static_cast<gcllPtrT>(nextChld)->getNoChld();

      switch (noChld)
      {
      case 0:
      {

         nextOk = nextChld->next;
         const size_t wasChild = getChildNo(nextChld, nextChld->parent);
         //std::cout << __LINE__ << " del  " << nextChld << "\n";
         //std::cout << __LINE__ << "      " << nextChld->parent << ":" << wasChild << " -> NULL \n";
         //std::cout << __LINE__ << " rew  " << curPtr << " -> " << nextOk << "\n";

         static_cast<gcllPtrT>(nextChld->parent)->child[wasChild] = NULL;
         delete static_cast<qcllPtrT>(nextChld);

         curPtr->next = nextOk;
         break;
      }

      case 1:
      {
         //
         // only delete chains not leading anywhere
         //
         std::cout << __LINE__ << " chai " << nextChld << " " << "\n";
         goNext();
         break;
         /*nextOk = nextChld->next;

         // follow the chain
         while (not nextOk->isParticle &&
                not nextOk->isCZ &&
                static_cast<gcllPtrT>(nextOk)->getNoChld() < 2)
         {
            nextOk = nextOk->next;
            std::cout << __LINE__ << " " << nextOk << "\n";
         }
         lastOk = curPtr;

         std::cout << __LINE__ << "\n";

         curPtr = lastOk->next;
         while (curPtr != nextOk)
         {
            const size_t wasChild = getChildNo(nextChld, nextChld->parent);
            static_cast<gcllPtrT>(nextChld->parent)->child[wasChild] = NULL;
            //std::cout << __LINE__ << " Xdelete " << nextChld << "\n";
            delete static_cast<qcllPtrT>(nextChld);
            //curPtr = curPtr->next;
            goNext();
         }
         std::cout << __LINE__ << "\n";

         lastOk->next = nextOk;

         break;*/
      }

      default:
      {
         std::cout << __LINE__ << " cell " << nextChld << " " << nextChld->isParticle << " " << nextChld->isCZ << " " << noChld << "\n";
         goNext();
         break;
      }
      }

      std::cout << __LINE__ << " lend " << curPtr << "\n";
   }

   std::cout << __LINE__ << " lastOk " << lastOk << " " << chldLastDummy << "\n";

   delete chldLastDummy;
   lastOk->next = chldLastNext;
   _czll->chldLast = lastOk;
   
   std::cout << __LINE__ << " chldLast " << chldLast << " " << curPtr << "\n";
   std::cout << __LINE__ << " finished\n\n\n";
}
};

#endif

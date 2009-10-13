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

#include "bhtree_worker.cpp"

namespace sphlatch {
class BHTreeHousekeeper : public BHTreeWorker {
public:

   BHTreeHousekeeper(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeHousekeeper(const BHTreeHousekeeper& _hk) : BHTreeWorker(_hk) { }
   ~BHTreeHousekeeper() { }

   ///
   /// set the next & skip pointers
   ///

   void setNext(const czllPtrT _czll);
   void setNextCZ();
   void setSkip();

private:
   void setNextRecursor();
   void setNextCZRecursor();

   nodePtrT lastPtr;

   std::vector<gcllPtrT> lastSkipeeAtDepth;
};

//FIXME: this is untested
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
   lastPtr->next = rootPtr;
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
   //FIXME: get maxdepth from tree
   //const size_t maxDepth = treePtr->maxDepth;
   const size_t maxDepth = 100;

   if (lastSkipeeAtDepth.size() != maxDepth)
      lastSkipeeAtDepth.resize(maxDepth);
   for (size_t i = 0; i < maxDepth; i++)
      lastSkipeeAtDepth[i] = NULL;

   ///
   /// do the next-walk
   ///
   while (curPtr != rootPtr)
   {
      goNext();
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
};

#endif

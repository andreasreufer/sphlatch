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
   void setPreorder(const czllPtrT _czll);


private:
   void setNextRecursor();

   nodePtrT lastPtr;
   size_t   lastDepth;

   std::vector<gcllPtrT> lastSkipeeAtDepth;
   size_t lastSkipeeDepth;
};


void BHTreeHousekeeper::setPreorder(const czllPtrT _czll)
{
   ///
   /// wire next pointer by doing a preorder tree walk
   ///
   curPtr  = _czll;
   lastPtr = _czll;
   setNextRecursor();
   lastPtr->next = rootPtr;

   curPtr = _czll;

   ///
   /// now do a preorder walk by using the next pointers
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
   //const size_t maxDepth = treePtr->maxDepth;
   const size_t maxDepth = 100;
   if ( lastSkipeeAtDepth.size() != maxDepth )
     lastSkipeeAtDepth.resize(maxDepth);
   for (size_t i = 0; i < maxDepth; i++)
     lastSkipeeAtDepth[i] = NULL;

   const nodePtrT stopPtr = _czll;

   lastSkipeeAtDepth[curPtr->depth] = _czll;
   lastSkipeeDepth = curPtr->depth;

   while( curPtr != stopPtr )
   {
     goNext();
     const size_t depth = curPtr->depth;
     if ( not curPtr->isParticle )
     {
       //lastSkipeeAtDepth[depth] 
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
         { }
      }
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
};

#endif

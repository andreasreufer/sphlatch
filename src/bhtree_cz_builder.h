#ifndef BHTREE_CZ_BUILDER_H
#define BHTREE_CZ_BUILDER_H

/*
 *  bhtree_cz_builder.h
 *
 *  Created by Andreas Reufer on 04.12.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_worker.h"
#include "communication_manager.h"
typedef sphlatch::CommunicationManager   commT;

using namespace sphlatch::vectindices;

namespace sphlatch {
class BHTreeCZBuilder : public BHTreeWorker {
public:

   //using BHTree::config;

   BHTreeCZBuilder(const treePtrT _treePtr) :
      BHTreeWorker(_treePtr),
      CommManager(commT::instance()),
      ///
      /// how to get this value cleverly?
      ///
      CZcosts(_treePtr->CZbottomCells.size())
   { }

   ~BHTreeCZBuilder()
   { }

   void operator()();

private:
   commT& CommManager;

   valvectType CZcosts;

   ///
   ///
   ///
   void refineCZcell(const czllPtrT _czllPtr);
   void gatherCZcell(const czllPtrT _czllPtr);
   void sumCZcosts();
};

///  - particles have been moved inside the tree
///  - do iterate CZ cells list
///    - globally sum up particle counts in CZ cells
///    - if a CZ cell contains too much particles, split it up
///      into CZ cell children and delete current CZ cell from list.
///      deleted normals cells which were replaced by CZ cells and rewire
///      everything again
///    - if a CZ cell contains too little particles, go to parent, add
///      up all its children and check if those are enough particles to
///      build a new CZ cell. if so, merge all children to the new CZ
///      cell and replace them by normal cells
///    while not all CZ cells are balanced
///
void BHTreeCZBuilder::operator()()
{
   bool CZbalanced = false;

   while (not CZbalanced)
   {
      czllPtrListItrT        CZlistItr = treePtr->CZbottomCells.begin();
      const czllPtrListCItrT CZlistEnd = treePtr->CZbottomCells.end();

      while (CZlistItr != CZlistEnd)
      {
///
///    - globally sum up particle counts in CZ cells
///    - if a CZ cell contains too much particles, split it up
///      into CZ cell children and delete current CZ cell from list.
///      delete normal cells which were replaced by CZ cells and rewire
///      everything again
///    - if a CZ cell contains too little particles, go to parent, add
///      up all its children and check if those are enough particles to
///      build a new CZ cell. if so, merge all children to the new CZ
///      cell and replace them by normal cells
///
         CZlistItr++;
      }
   }

   czllPtrListCItrT       CZlistItr = treePtr->CZbottomCells.begin();
   const czllPtrListCItrT CZlistEnd = treePtr->CZbottomCells.end();

   while (CZlistItr != CZlistEnd)
   {
      curPtr = (*CZlistItr);
      CZlistItr++;
   }
}

void BHTreeCZBuilder::sumCZcosts()
{
   ///
   /// CZ cell list -> vector
   ///
   czllPtrListCItrT       CZlistCItr = treePtr->CZbottomCells.begin();
   const czllPtrListCItrT CZlistEnd  = treePtr->CZbottomCells.end();
   size_t i = 0;

   while (CZlistCItr != CZlistEnd)
   {
      CZcosts[i] = (*CZlistCItr)->cost;
      i++;
   }

   ///
   /// sum up vector
   ///
   CommManager.sum(CZcosts);

   ///
   /// vector -> CZ cell list
   ///
   czllPtrListItrT CZlistItr = treePtr->CZbottomCells.begin();
   i = 0;

   while (CZlistItr != CZlistEnd)
   {
      (*CZlistItr)->cost = CZcosts[i];
      i++;
   }
}

void BHTreeCZBuilder::refineCZcell(const czllPtrT _czllPtr)
{
   ///
   /// transform all children to empty CZ cells. if a child
   /// is a particle, first link a empty cell in between.
   ///
   curPtr = _czllPtr;
   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
      {
         goChild(i);
         if (curPtr->isParticle)
         {
           const partPtrT resPartPtr = static_cast<partPtrT>(curPtr);
           curPtr = treePtr->cellAllocator.pop();
          static_cast<cellPtrT>(curPtr)->clear();
          static_cast<cellPtrT>(curPtr)->inheritCellPos(i);
          curPtr->parent = _czllPtr;

          resPartPtr->parent = curPtr;
          //pushDownSingle(resPartPtr);

            /*resPartPtr = static_cast<gcllPtr>(curPtr)->child[i];
            static_cast<gcllPtr>(curPtr)->child[i] = 
              t*/
         }

         //static_cast<gcllPtr>(curPtr)->

         /*if (
            if ( static_cast<gcllPtr>(curPtr)->child[i] != NULL )
            {
            }
            //if ( _czllPtr->child[i]->isParticle )
            if ( static_cast<gcllPtr>(curPtr)->child[i]->isParticle )
            {
            const partPtrT resPartPtr = static_cast<partPtrT>(_czllPtr->child[i]);
            _czllPtr->child[i] = treePtr->cellAllocator.pop();

            }*/
      }
   }

   treePtr->CZbottomCells.erase(_czllPtr->listItr);
   treePtr->CZbottomCells.push_front(_czllPtr);
}

void BHTreeCZBuilder::gatherCZcell(const czllPtrT _czllPtr)
{
   treePtr->CZbottomCells.erase(_czllPtr->listItr);
   treePtr->CZbottomCells.push_front(_czllPtr);
}
}

#endif

#ifndef BHTREE_CZ_BUILDER_CPP
#define BHTREE_CZ_BUILDER_CPP

/*
 *  bhtree_cz_builder.cpp
 *
 *  Created by Andreas Reufer on 04.12.08.
 *  Copyright 2008 University of Berne. All rights reserved.
 *
 */

#include "typedefs.h"
#include "bhtree_worker.cpp"
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
      CZcosts(_treePtr->maxCZBottCells)
   { }

   ~BHTreeCZBuilder()
   { }

   void operator()();

private:
   commT& CommManager;

   fvectT CZcosts;

   ///
   ///
   ///
   void refineCZcell(const czllPtrT _czllPtr);
   void gatherCZcell(const czllPtrT _czllPtr);
   void sumCZcosts();
   void sumCostRecursor();
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

   const fType costHighMark = 1000.;
   const fType costLowMark  = 10.;

   std::cout << "start to balance ...\n";
   while (not CZbalanced)
   {
      czllPtrListT::iterator       CZlistItr = treePtr->CZbottom.begin();
      czllPtrListT::const_iterator CZlistEnd = treePtr->CZbottom.end();

      while (CZlistItr != CZlistEnd)
      {
         ///
         /// sum up CZ cost
         ///
         sumCZcosts();

         ///
         ///    - if a CZ cell contains too much particles, split it up
         ///      into CZ cell children and delete current CZ cell from list.
         ///      delete normal cells which were replaced by CZ cells and rewire
         ///      everything again
         ///    - if a CZ cell contains too little particles, go to parent, add
         ///      up all its children and check if those are enough particles to
         ///      build a new CZ cell. if so, merge all children to the new CZ
         ///      cell and replace them by normal cells
         ///

         if ((*CZlistItr)->absCost > costHighMark)
         {
            refineCZcell(*CZlistItr);
            CZbalanced = false;
         }
         else
         if ((*CZlistItr)->parent != NULL)
         {
            if (static_cast<czllPtrT>((*CZlistItr)->parent)->absCost <
                costLowMark)
               gatherCZcell(static_cast<czllPtrT>((*CZlistItr)->parent));
            CZbalanced = false;
         }
         CZlistItr++;
      }
      CZbalanced = true;
   }

   /*czllPtrListCItrT       CZlistItr = treePtr->CZbottomCells.begin();
      const czllPtrListCItrT CZlistEnd = treePtr->CZbottomCells.end();

      while (CZlistItr != CZlistEnd)
      {
      curPtr = (*CZlistItr);
      CZlistItr++;
      }*/
}

void BHTreeCZBuilder::sumCZcosts()
{
   ///
   /// CZ cell list -> vector
   ///
   czllPtrListT::const_iterator CZlistCItr = treePtr->CZbottom.begin();
   czllPtrListT::const_iterator CZlistEnd  = treePtr->CZbottom.end();
   size_t i = 0;

   while (CZlistCItr != CZlistEnd)
   {
      CZcosts[i] = (*CZlistCItr)->absCost;
      i++;
      CZlistCItr++;
   }

   ///
   /// sum up vector
   ///
   CommManager.sum(CZcosts);

   ///
   /// vector -> CZ cell list
   ///
   czllPtrListT::iterator CZlistItr = treePtr->CZbottom.begin();
   i = 0;

   while (CZlistItr != CZlistEnd)
   {
      (*CZlistItr)->absCost = CZcosts[i];
      i++;
      CZlistItr++;
   }

   ///
   /// make sure the costs of all CZ cells is added up
   ///
std::cout << __LINE__ << "\n";
   goRoot();
std::cout << __LINE__ << "\n";
   sumCostRecursor();
std::cout << __LINE__ << "\n";
}

//FIXME: untested
void BHTreeCZBuilder::sumCostRecursor()
{
   if (static_cast<gcllPtrT>(curPtr)->atBottom)
   {
std::cout << __LINE__ << "\n";
      if (curPtr->parent != NULL)
         static_cast<czllPtrT>(curPtr->parent)->absCost +=
            static_cast<czllPtrT>(curPtr)->absCost;
std::cout << __LINE__ << "\n";
   }
   else
   {
std::cout << __LINE__ << "\n";
      static_cast<czllPtrT>(curPtr)->absCost = 0.;
std::cout << __LINE__ << "\n";
      for (size_t i = 0; i < 8; i++)
      {
std::cout << __LINE__ << "\n";
         goChild(i);
std::cout << __LINE__ << "\n";
         sumCostRecursor();
std::cout << __LINE__ << "\n";
         goUp();
std::cout << __LINE__ << "\n";
      }
std::cout << __LINE__ << "\n";
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
            const pnodPtrT resPartPtr = static_cast<pnodPtrT>(curPtr);
            //curPtr = treePtr->cellAllocator.pop();
            curPtr = new czllT;
            static_cast<gcllPtrT>(curPtr)->clear();
            static_cast<gcllPtrT>(curPtr)->inheritCellPos(i);
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

   /*treePtr->CZbottomCells.erase(_czllPtr->listItr);
      treePtr->CZbottomCells.push_front(_czllPtr);*/
}

void BHTreeCZBuilder::gatherCZcell(const czllPtrT _czllPtr)
{
   /*treePtr->CZbottomCells.erase(_czllPtr->listItr);
      treePtr->CZbottomCells.push_front(_czllPtr);*/
}
};

#endif

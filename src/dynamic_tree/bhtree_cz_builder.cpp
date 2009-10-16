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

#include "bhtree_treedump.cpp"
typedef sphlatch::BHTreeDump               dumpT;

using namespace sphlatch::vectindices;

namespace sphlatch {
class BHTreeCZBuilder : public BHTreeWorker {
public:
   BHTreeCZBuilder(const treePtrT _treePtr) :
      BHTreeWorker(_treePtr),
      CommManager(commT::instance()),
      CZcosts(_treePtr->maxCZBottCells),
      CZparts(_treePtr->maxCZBottCells),
      dumper(_treePtr)
   { }

   ~BHTreeCZBuilder()
   { }

   void operator()();

private:
   commT& CommManager;

   fvectT CZcosts;
   ivectT CZparts;


   ///
   ///
   ///
   void refineCZcell(const czllPtrT _czllPtr);

   void gatherCZcell(const czllPtrT _czllPtr);

   void sumCZcosts();
   void sumCostRecursor();

   void countPartsRecursor();

   size_t countParts;

   dumpT dumper;
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
   std::list<czllPtrT>& CZbottom(treePtr->CZbottom);

   bool CZbalanced = false;

   const fType costHighMark = 5.;
   const fType costLowMark  = 2.;

   std::cout << "start to balance ...\n";
   while (not CZbalanced)
   {
      ///
      /// sum up CZ cost
      ///
      sumCZcosts();

      std::cout << __LINE__ << "\n";
      czllPtrListT::iterator       CZlistItr = CZbottom.begin();
      czllPtrListT::const_iterator CZlistEnd = CZbottom.end();
      std::cout << __LINE__ << "\n";

      //while (CZlistItr != CZlistEnd)
      {
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

         CZbalanced = true;
         if ((*CZlistItr)->absCost > costHighMark)
         {
            const czllPtrT refdCellPtr = *CZlistItr;
            // refine the cell
            std::cout << __LINE__ << " refine cell " << refdCellPtr << " " << refdCellPtr->atBottom << "\n";
            refineCZcell(refdCellPtr);
            std::cout << __LINE__ << " refined cell\n";

            // delete from list and insert new ones
            for (size_t i = 0; i < 8; i++)
            {
               if (refdCellPtr->child[i] != NULL)
               {
                  std::cout << __LINE__ << ":" << i << "\n";
                  CZbottom.insert(CZlistItr,
                                  static_cast<czllPtrT>(refdCellPtr->child[i]));
                  static_cast<czllPtrT>(refdCellPtr->child[i])->atBottom = true;
               }
            }

            refdCellPtr->atBottom = false;
            CZbottom.erase(CZlistItr);

            CZbalanced = false;
         }
         else
         if ((*CZlistItr)->parent != NULL)
         {
            goUp();
            if (static_cast<czllPtrT>(curPtr)->absCost < costLowMark)
               gatherCZcell(static_cast<czllPtrT>(curPtr));
            CZbalanced = false;
         }
         std::cout << __LINE__ << " next CZlist \n";
         CZlistItr++;
      }
   }
   std::cout << __LINE__ << "\n";
}

void BHTreeCZBuilder::sumCZcosts()
{
   std::cout << __LINE__ << "sumCZcosts()\n";
   ///
   /// CZ cell list -> vector
   ///
   czllPtrListT::const_iterator CZlistCItr = treePtr->CZbottom.begin();
   czllPtrListT::const_iterator CZlistEnd  = treePtr->CZbottom.end();
   size_t i = 0;

   while (CZlistCItr != CZlistEnd)
   {
      CZcosts[i] = (*CZlistCItr)->absCost;
      CZparts[i] = (*CZlistCItr)->noParts;
      i++;
      CZlistCItr++;
   }

   ///
   /// sum up vector
   ///
   CommManager.sum(CZcosts);
   CommManager.sumUpCounts(CZparts);

   ///
   /// vector -> CZ cell list
   ///
   czllPtrListT::iterator CZlistItr = treePtr->CZbottom.begin();
   i = 0;

   while (CZlistItr != CZlistEnd)
   {
      (*CZlistItr)->absCost = CZcosts[i];
      (*CZlistItr)->noParts = CZparts[i];
      i++;
      CZlistItr++;
   }

   ///
   /// make sure the costs of all CZ cells is added up
   ///
   std::cout << __LINE__ << "\n";
   goRoot();
   sumCostRecursor();
   std::cout << __LINE__ << "\n";
}

//FIXME: untested
void BHTreeCZBuilder::sumCostRecursor()
{
   std::cout << __LINE__ << ":" << curPtr << "\n";
   if (static_cast<gcllPtrT>(curPtr)->atBottom)
   {
      if (curPtr->parent != NULL)
      {
         static_cast<czllPtrT>(curPtr->parent)->absCost +=
            static_cast<czllPtrT>(curPtr)->absCost;
         static_cast<czllPtrT>(curPtr->parent)->noParts +=
            static_cast<czllPtrT>(curPtr)->noParts;
      }
   }
   else
   {
      static_cast<czllPtrT>(curPtr)->absCost = 0.;
      static_cast<czllPtrT>(curPtr)->noParts = 0;
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            std::cout << __LINE__ << ":" << curPtr->depth << "\n";
            goChild(i);
            sumCostRecursor();
            goUp();
         }
      }
   }
}

void BHTreeCZBuilder::refineCZcell(const czllPtrT _czllPtr)
{
   std::cout << "refine cell!\n";
   std::cout << __LINE__ << " refine\n";
   ///
   /// transform all children to empty CZ cells. if a child
   /// is a particle, first link a empty cell in between.
   ///
   curPtr = _czllPtr;

   const fType relCost = (_czllPtr->absCost / _czllPtr->noParts);
   std::cout << __LINE__ << "\n";

   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
      {
         goChild(i);
         std::cout << __LINE__ << "\n";

         if (curPtr->isParticle)
         {
            std::cout << __LINE__ << "\n";
            const pnodPtrT resPartPtr = static_cast<pnodPtrT>(curPtr);

            curPtr = new czllT;
            static_cast<gcllPtrT>(curPtr)->clear();
            static_cast<gcllPtrT>(curPtr)->inheritCellPos(i);
            curPtr->parent = _czllPtr;

            const size_t newOct = getOctant(resPartPtr->pos);
            static_cast<gcllPtrT>(curPtr)->child[newOct] = resPartPtr;
            resPartPtr->depth++;
            resPartPtr->parent = curPtr;
         }
         else if (not curPtr->isCZ)
         {
            std::cout << __LINE__ << "\n";
            const czllPtrT newCZllPtr = new czllT;
            newCZllPtr->initFromCell(*static_cast<czllPtrT>(curPtr));

            delete static_cast<qcllPtrT>(curPtr);
            curPtr = newCZllPtr;
         }

         ///
         /// now count the particles of this child
         /// estimate the absolute cost of those
         ///
         std::cout << __LINE__ << "\n";
         countParts = 0;
         countPartsRecursor();
         static_cast<czllPtrT>(curPtr)->absCost = countParts * relCost;
         static_cast<czllPtrT>(curPtr)->noParts = countParts;
         std::cout << __LINE__ << "\n";

         goUp();
         std::cout << __LINE__ << "\n";
      }
   }
   std::cout << __LINE__ << "refine finished\n";
}

void BHTreeCZBuilder::countPartsRecursor()
{
   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i]->isParticle)
            countParts++;
         else
         {
            goChild(i);
            countPartsRecursor();
            goUp();
         }
      }
   }
}

void BHTreeCZBuilder::gatherCZcell(const czllPtrT _czllPtr)
{
   std::cout << "gather!\n";
   for (size_t i = 0; i < 8; i++)
   {
      if (_czllPtr->child[i] != NULL)
      {
         if (_czllPtr->child[i]->atBottom)
         {
            const qcllPtrT newCell = new qcllT;
            const czllPtrT oldCell = static_cast<czllPtrT>(_czllPtr->child[i]);

            newCell->initFromCZll(*oldCell);
            treePtr->CZbottom.remove(oldCell);
            delete oldCell;
         }
         else
            gatherCZcell(static_cast<czllPtrT>(_czllPtr->child[i]));
      }
   }
}
};

#endif

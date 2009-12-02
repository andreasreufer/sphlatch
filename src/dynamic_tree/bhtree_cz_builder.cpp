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

#ifdef SPHLATCH_MPI
#include "communication_manager.h"
typedef sphlatch::CommunicationManager   commT;
#endif

#include "bhtree_treedump.cpp"
typedef sphlatch::BHTreeDump             dumpT;

namespace sphlatch {
class BHTreeCZBuilder : public BHTreeWorker {
public:
   BHTreeCZBuilder(const treePtrT _treePtr) :
      BHTreeWorker(_treePtr),
#ifdef SPHLATCH_MPI
      CommManager(commT::instance()),
#endif
      CZcosts(_treePtr->maxCZBottCells),
      CZparts(_treePtr->maxCZBottCells),
      gathOrphFrst(NULL),
      gathOrphLast(NULL),
      dumper(_treePtr)
   { }

   ~BHTreeCZBuilder()
   { }

   void rebalance(const fType _lowMark, const fType _highMark);

private:
#ifdef SPHLATCH_MPI
   commT& CommManager;
#endif

   fvectT CZcosts;
   ivectT CZparts;

   void refineCZcell(const czllPtrT _czllPtr);
   void gatherCZcell(const czllPtrT _czllPtr);

   nodePtrT gathOrphFrst, gathOrphLast;

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
void BHTreeCZBuilder::rebalance(const fType _lowMark, const fType _highMark)
{
   std::list<czllPtrT>& CZbottom(treePtr->CZbottom);

   bool CZbalanced = false;

   const fType costHighMark = _highMark;
   const fType costLowMark  = _lowMark;

   while (not CZbalanced)
   {
      ///
      /// sum up CZ cost
      ///
      sumCZcosts();

      czllPtrListT::iterator       CZlistItr = CZbottom.begin();
      czllPtrListT::const_iterator CZlistEnd = CZbottom.end();

      CZbalanced = true;
      while (CZlistItr != CZlistEnd)
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

         if ((*CZlistItr)->absCost > costHighMark)
         {
            CZbalanced = false;
            const czllPtrT refdCellPtr = *CZlistItr;

            // refine the cell
            refineCZcell(refdCellPtr);

            // delete from list and insert new ones
            for (int i = 7; i >= 0; i--)
            {
               CZbottom.insert(CZlistItr,
                               static_cast<czllPtrT>(refdCellPtr->child[i]));
               static_cast<czllPtrT>(refdCellPtr->child[i])->atBottom = true;
            }
            refdCellPtr->atBottom = false;

            CZlistItr--;
            CZbottom.remove(refdCellPtr);
         }
         else
         if ((*CZlistItr)->parent != NULL)
         {
            const czllPtrT gathCellPtr =
               static_cast<czllPtrT>((*CZlistItr)->parent);
            if (gathCellPtr->absCost < costLowMark)
            {
               ///
               /// go back again in the CZ cell list just to make sure,
               /// that the iterator doesn't point to a CZ cell which was
               /// deleted by gatherCZcell()
               ///
               CZlistItr--;

               gathOrphFrst = NULL;
               gathOrphLast = NULL;
               gatherCZcell(gathCellPtr);
               gathCellPtr->orphFrst = static_cast<pnodPtrT>(gathOrphFrst);
               gathCellPtr->orphLast = static_cast<pnodPtrT>(gathOrphLast);

               gathCellPtr->atBottom = true;
               CZbottom.insert(CZlistItr, gathCellPtr);

               CZbalanced = false;
            }
         }
         CZlistItr++;
      }
   }
}

// FIXME: move this to tree
void BHTreeCZBuilder::sumCZcosts()
{
#ifdef SPHLATCH_MPI
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
#endif
   ///
   /// make sure the costs of all CZ cells is added up
   ///
   goRoot();
   sumCostRecursor();

   goRoot();
}

void BHTreeCZBuilder::sumCostRecursor()
{
   if (not curPtr->atBottom)
   {
      static_cast<czllPtrT>(curPtr)->absCost = 0.;
      static_cast<czllPtrT>(curPtr)->noParts = 0;

      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            goChild(i);
            sumCostRecursor();
            static_cast<czllPtrT>(curPtr->parent)->absCost +=
               static_cast<czllPtrT>(curPtr)->absCost;
            static_cast<czllPtrT>(curPtr->parent)->noParts +=
               static_cast<czllPtrT>(curPtr)->noParts;
            goUp();
         }
      }
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
      if (static_cast<gcllPtrT>(curPtr)->child[i] == NULL)
      {
         static_cast<gcllPtrT>(curPtr)->child[i] = new czllT;
         goChild(i);
         static_cast<czllPtrT>(curPtr)->clear();
         curPtr->parent = _czllPtr;
         static_cast<gcllPtrT>(curPtr)->inheritCellPos(i);
         goUp();
      }
      else if (static_cast<gcllPtrT>(curPtr)->child[i]->isParticle)
      {
         const pnodPtrT resPartPtr =
            static_cast<pnodPtrT>(static_cast<gcllPtrT>(curPtr)->child[i]);

         static_cast<gcllPtrT>(curPtr)->child[i] = new czllT;

         goChild(i);
         static_cast<czllPtrT>(curPtr)->clear();
         curPtr->parent = _czllPtr;
         static_cast<czllPtrT>(curPtr)->inheritCellPos(i);

         const size_t newOct = getOctant(resPartPtr->pos);
         static_cast<gcllPtrT>(curPtr)->child[newOct] = resPartPtr;
         resPartPtr->depth++;
         resPartPtr->parent = curPtr;
         goUp();
      }
      else
      if (not static_cast<gcllPtrT>(curPtr)->child[i]->isCZ)
      {
         const qcllPtrT oldCell =
            static_cast<qcllPtrT>(static_cast<gcllPtrT>(curPtr)->child[i]);
         static_cast<gcllPtrT>(curPtr)->child[i] = new czllT;
         goChild(i);
         static_cast<czllPtrT>(curPtr)->clear();
         static_cast<czllPtrT>(curPtr)->initFromCell(*oldCell);
         delete oldCell;
         goUp();
      }
   }


   ///
   /// now count the particles of this child
   /// estimate the absolute cost of those
   ///
   const fType relCost = (_czllPtr->absCost / _czllPtr->noParts);
   for (size_t i = 0; i < 8; i++)
   {
      goChild(i);
      countParts = 0;
      countPartsRecursor();
      static_cast<czllPtrT>(curPtr)->absCost = countParts * relCost;
      static_cast<czllPtrT>(curPtr)->noParts = countParts;

      goUp();
   }

   ///
   /// now also distribute the orphans of the former parent cell to
   /// the new bottom cells
   ///
   nodePtrT curOrph = _czllPtr->orphFrst;
   nodePtrT nxtOrph = NULL;
   curPtr = _czllPtr;

   while (curOrph != NULL)
   {
      // next ptr will be overwritten in adopt(), so we need a temporary
      // storage for the next orphan
      nxtOrph = curOrph->next;

      const size_t newOct = getOctant(static_cast<pnodPtrT>(curOrph)->pos);
      static_cast<czllPtrT>(_czllPtr->child[newOct])->adopt(
         static_cast<pnodPtrT>(curOrph));
      static_cast<czllPtrT>(_czllPtr->child[newOct])->noParts++;
      static_cast<czllPtrT>(_czllPtr->child[newOct])->absCost += relCost;

      curOrph = nxtOrph;
   }

   static_cast<czllPtrT>(curPtr)->orphFrst = NULL;
   static_cast<czllPtrT>(curPtr)->orphLast = NULL;
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

//FIXME: orphans are not collected!
void BHTreeCZBuilder::gatherCZcell(const czllPtrT _czllPtr)
{
   for (size_t i = 0; i < 8; i++)
   {
      if (_czllPtr->child[i] != NULL)
      {
         if (_czllPtr->child[i]->atBottom)
         {
            const czllPtrT oldCell = static_cast<czllPtrT>(_czllPtr->child[i]);

            if (oldCell->orphFrst != NULL)
            {
               if (gathOrphFrst == NULL)
               {
                  gathOrphFrst = oldCell->orphFrst;
                  gathOrphLast = oldCell->orphLast;
               }
               else
               {
                  gathOrphLast->next = oldCell->orphFrst;
                  gathOrphLast       = oldCell->orphLast;
               }
            }


            treePtr->CZbottom.remove(oldCell);
            _czllPtr->child[i] = czllToCell(oldCell);
         }
         else
         {
            gatherCZcell(static_cast<czllPtrT>(_czllPtr->child[i]));
         }
      }
   }
}
};

#endif

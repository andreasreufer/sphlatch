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
typedef sphlatch::BHTreeDump             dumpT;

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

   void rebalance(const fType _lowMark, const fType _highMark);

private:
   commT& CommManager;

   fvectT CZcosts;
   ivectT CZparts;

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
void BHTreeCZBuilder::rebalance(const fType _lowMark, const fType _highMark)
{
   std::list<czllPtrT>& CZbottom(treePtr->CZbottom);

   bool CZbalanced = false;

   const fType costHighMark = _highMark;
   const fType costLowMark  = _lowMark;

   size_t round = 1;
   while (not CZbalanced)
   {
      std::ostringstream roundStr;
      roundStr << round;
      std::string dumpName = "dumpCZ";
      dumpName.append(roundStr.str());

      dumper.dotDump(dumpName + ".dot");
      dumper.ptrDump(dumpName + ".txt");
      //std::cout << "ROUND " << round << "\n";
      round++;

      ///
      /// sum up CZ cost
      ///
      sumCZcosts();
      //std::cout << "summed up CZ\n\n";

      czllPtrListT::iterator       CZlistItr = CZbottom.begin();
      czllPtrListT::const_iterator CZlistEnd = CZbottom.end();

      while (CZlistItr != CZlistEnd)
      {
        const czllPtrT curCZ = *CZlistItr;
        //std::cout << curCZ << " " << curCZ->absCost << " " << curCZ->depth << "\n";
        CZlistItr++;
      }
      //std::cout << "\n";

      CZlistItr = CZbottom.begin();

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

         CZbalanced = true;
         if ((*CZlistItr)->absCost > costHighMark)
         {
            CZbalanced = false;
            const czllPtrT refdCellPtr = *CZlistItr;
            
            // refine the cell
            //std::cout << "refine " << refdCellPtr << " " << refdCellPtr->absCost << "\n";
            refineCZcell(refdCellPtr);

            // delete from list and insert new ones
            for (int i = 7; i >= 0; i--)
            {
               CZbottom.insert(CZlistItr,
                               static_cast<czllPtrT>(refdCellPtr->child[i]));
               static_cast<czllPtrT>(refdCellPtr->child[i])->atBottom = true;
               //std::cout << "add new " << refdCellPtr->child[i] << "\n";
            }
            refdCellPtr->atBottom = false;

            CZlistItr--;
            CZbottom.remove(refdCellPtr);
            //std::cout << "rem old " << refdCellPtr << "\n";
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
               //std::cout << "gath call " << gathCellPtr << "\n";
               gatherCZcell(gathCellPtr);

               gathCellPtr->atBottom = true;
               CZbottom.insert(CZlistItr, gathCellPtr);

               CZbalanced = false;
            }
         }
         CZlistItr++;
         //std::cout << __LINE__ << " " << *CZlistItr - *CZlistFrst << "\n";
      }
   }
   //std::cout << __LINE__ << "\n";
}

// FIXME: move this to tree
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
   
   //std::cout << "refine " << _czllPtr << "   " << _czllPtr->orphFrst << " " << _czllPtr->orphLast << " " << _czllPtr->orphFrst << "\n";

   //curPtr = _czllPtr->orphFrst;

   pnodPtrT curOrph = _czllPtr->orphFrst;
   curPtr = _czllPtr;
   while (curOrph != NULL)
   {
     const size_t newOct = getOctant(curOrph->pos);

     static_cast<czllPtrT>(_czllPtr->child[newOct])->adopt(curOrph);
      static_cast<czllPtrT>(_czllPtr->child[newOct])->noParts++;
      static_cast<czllPtrT>(_czllPtr->child[newOct])->absCost += relCost;

     std::cout << "next orphan " << static_cast<pnodPtrT>(curOrph->next);
     curOrph = static_cast<pnodPtrT>(curOrph->next);
   }

   /*while (curPtr != NULL)
   {
      const pnodPtrT curPart = static_cast<pnodPtrT>(curPtr);
      const size_t   newOct  = getOctant(curPart->pos);
      std::cout << curPtr << " (" << curPtr->ident << ") move to " 
                << static_cast<czllPtrT>(_czllPtr->child[newOct])
                << "\n";
      static_cast<czllPtrT>(_czllPtr->child[newOct])->adopt(curPart);
      static_cast<czllPtrT>(_czllPtr->child[newOct])->noParts++;
      static_cast<czllPtrT>(_czllPtr->child[newOct])->absCost += relCost;
      goNext();
   }
   curPtr = _czllPtr;*/
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

            std::cout << "rem. oldCell " << oldCell->orphFrst << " " << oldCell->orphLast << "\n";

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

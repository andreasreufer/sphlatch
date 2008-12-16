#ifndef BHTREE_PART_INSERTER_H
#define BHTREE_PART_INSERTER_H

/*
 *  bhtree_part_inserter.h
 *
 *  Created by Andreas Reufer on 02.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"
#include "bhtree_errhandler_tmp.h"

namespace sphlatch {
class BHTreePartsInserter : public BHTreeWorker {
public:

   typedef BHTreePartsInserter   selfT;

   BHTreePartsInserter(treePtrT _treePtr) :  BHTreeWorker(_treePtr) { }
   BHTreePartsInserter(const selfT& _inserter) : BHTreeWorker(_inserter) { }
   ~BHTreePartsInserter() { }

public:
   void operator()(const size_t _i);

private:
   void recursor(const size_t _i);
};

///
/// entry function to insert particle:
/// - check whether particle lies inside root cell
/// - if so, call insertion recursor
///
void BHTreePartsInserter::operator()(const size_t _i)
{
   goRoot();

   ///
   /// check whether the particles comes to lie in the root
   /// cell, if not throw an exception
   ///
   const fType curX = pos(_i, X);
   const fType curY = pos(_i, Y);
   const fType curZ = pos(_i, Z);

   if (not pointInsideCell(curX, curY, curZ))
      throw PartOutsideTree(_i, rootPtr);

   recursor(_i);
}

///
/// particle insertion recursor
///
///
void BHTreePartsInserter::recursor(const size_t _i)
{
   const fType curX    = pos(_i, X);
   const fType curY    = pos(_i, Y);
   const fType curZ    = pos(_i, Z);
   const fType curCost = cost(_i);

   const size_t tarOct = getOctant(curX, curY, curZ);


   ///
   /// If targeted child is empty, place the particle there
   ///
   if (static_cast<cellPtrT>(curPtr)->child[tarOct] == NULL)
   {
      goChild(tarOct);

      /*curNodePtr->isEmpty = false;

         newPartChild(targetOctant);
         goChild(targetOctant);

         if (_newPartIdx >= noLocParts)
         {
          curNodePtr->isLocal = false;
         }
         else
         {
          curNodePtr->isLocal = true;
          /// save the particle's proxy
          partProxies[_newPartIdx] = curNodePtr;
         }

         /// particle saves its position to node directly
         static_cast<partPtrT>(curNodePtr)->xPos = pos(_newPartIdx, X);
         static_cast<partPtrT>(curNodePtr)->yPos = pos(_newPartIdx, Y);
         static_cast<partPtrT>(curNodePtr)->zPos = pos(_newPartIdx, Z);
         static_cast<partPtrT>(curNodePtr)->mass = m(_newPartIdx);

         /// ident saves the rowIndex of the particle
         curNodePtr->ident = _newPartIdx;

         goUp();*/
   }

   ///
   /// ... or if existing child is a node, then try to place the particle
   /// as a child of this node
   ///
   else if (not static_cast<cellPtrT>(curPtr)->child[tarOct]->isParticle)
   {
      goChild(tarOct);

      /*
         curNodePtr->isEmpty = false;
         goChild(targetOctant);
         insertParticleRecursor(_newPartIdx);
         goUp();*/
   }

   ///
   /// ... or if the existing child is a particle, then replace it by a
   /// cell node and try to insert them both again at the current node
   ///
   else
   {
      /*
         /// goto child, save resident particle
         goChild(targetOctant);
         size_t residentIdx = static_cast<partPtrT>(curNodePtr)->ident;
         goUp();

         /// replace particle by cell node
         delete static_cast<cellPtrT>(curNodePtr)->child[targetOctant];
         static_cast<cellPtrT>(curNodePtr)->child[targetOctant] = NULL;
         newCellChild(targetOctant);

         /// and try to insert both particles again
         goChild(targetOctant);
         ///
         /// check whether we are too deep
         ///
         if (curNodePtr->depth > 128)
         throw PartsTooClose(curNodePtr->depth,
                            residentIdx,
                            _newPartIdx,
                            rootPtr);

         insertParticleRecursor(residentIdx);
         insertParticleRecursor(_newPartIdx);
         goUp();*/
   }
}
};

#endif

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
  
   BHTreeCZBuilder(treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeCZBuilder(const BHTreeCZBuilder& _dumper) : BHTreeWorker(_dumper) { }

   /*BHTreeCZBuilder(BHTree* _treePtr) :
     CommManager(commT::instance())
   {
      /// get this from the tree
      CZpCounts.resize(4096);
      treePtr = _treePtr;
   }*/

   ~BHTreeCZBuilder()
   {
   }

   void operator()();

private:
   czllPtrlistT::iterator CZlistItr, CZlistEnd;
   countsVectType CZpCounts;

   commT& CommManager;

   void sumCZpCounts();
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
      CZlistItr = treePtr->CZbottomCells.begin();
      CZlistEnd = treePtr->CZbottomCells.end();

      while (CZlistItr != CZlistEnd)
      {
///
///    - globally sum up particle counts in CZ cells
///    - if a CZ cell contains too much particles, split it up
///      into CZ cell children and delete current CZ cell from list.
///      deleted normals cells which were replaced by CZ cells and rewire
///      everything again
///    - if a CZ cell contains too little particles, go to parent, add
///      up all its children and check if those are enough particles to
///      build a new CZ cell. if so, merge all children to the new CZ
///      cell and replace them by normal cells
///
         CZlistItr++;
      }
   }

   CZlistItr = treePtr->CZbottomCells.begin();
   CZlistEnd = treePtr->CZbottomCells.end();

   while (CZlistItr != CZlistEnd)
   {
      curPtr = (*CZlistItr);
      CZlistItr++;
   }
}

void BHTreeCZBuilder::sumCZpCounts()
{
   ///
   /// CZ cell list -> vector
   ///
   CZlistItr = treePtr->CZbottomCells.begin();
   CZlistEnd = treePtr->CZbottomCells.end();
   size_t i = 0;

   while (CZlistItr != CZlistEnd)
   {
      CZpCounts[i] = (*CZlistItr)->nParts;
      i++;
   }

   ///
   /// sum up vector
   ///
   CommManager.sumUpCounts(CZpCounts);

   ///
   /// vector -> CZ cell list
   ///
   CZlistItr = treePtr->CZbottomCells.begin();
   CZlistEnd = treePtr->CZbottomCells.end();
   i         = 0;

   while (CZlistItr != CZlistEnd)
   {
      (*CZlistItr)->nParts = CZpCounts[i];
      i++;
   }
}
};

#endif

/*
 * public:
 * ///
 * /// method to insert particle:
 * /// - go to root
 * /// - call the insertion recursor
 * ///
 * void insertParticle(size_t _newPartIdx)
 * {
 * ///
 * /// check whether the particle comes to lie in the root cell
 * ///
 * const fType curX = pos(_newPartIdx, X);
 * const fType curY = pos(_newPartIdx, Y);
 * const fType curZ = pos(_newPartIdx, Z);
 * const fType rootSize = static_cast<cellPtrT>(rootPtr)->cellSize;
 * const fType rootX = static_cast<cellPtrT>(rootPtr)->xCenter;
 * const fType rootY = static_cast<cellPtrT>(rootPtr)->yCenter;
 * const fType rootZ = static_cast<cellPtrT>(rootPtr)->zCenter;
 *
 * if (curX < rootX - 0.5 * rootSize ||
 *    curX > rootX + 0.5 * rootSize ||
 *    curY < rootY - 0.5 * rootSize ||
 *    curY > rootY + 0.5 * rootSize ||
 *    curZ < rootZ - 0.5 * rootSize ||
 *    curZ > rootZ + 0.5 * rootSize)
 *  throw PartOutsideTree(_newPartIdx, rootPtr);
 *
 * ///
 * /// everything seems to be fine, so let's try
 * /// to insert the particle
 * ///
 * goRoot();
 * insertParticleRecursor(_newPartIdx);
 *
 * if (_newPartIdx >= noLocParts)
 *  {
 *    ghostCounter++;
 *  }
 * else
 *  {
 *    partCounter++;
 *  }
 * }
 *
 * private:
 * ///
 * /// recursor for inserting a new particle:
 * /// try to insert as child of current
 * /// node. if child is
 * /// - empty, insert particle. we're done.
 * /// - a node, go to child and call recursor
 * /// - a particle, disconnect particle and
 * ///   call recursor for existing and new
 * ///   particle.
 * ///
 * void insertParticleRecursor(const size_t _newPartIdx)
 * {
 * size_t targetOctant = 0;
 *
 * targetOctant += (pos(_newPartIdx, X) <
 *                 static_cast<cellPtrT>(curNodePtr)->xCenter) ? 0 : 1;
 * targetOctant += (pos(_newPartIdx, Y) <
 *                 static_cast<cellPtrT>(curNodePtr)->yCenter) ? 0 : 2;
 * targetOctant += (pos(_newPartIdx, Z) <
 *                 static_cast<cellPtrT>(curNodePtr)->zCenter) ? 0 : 4;
 *
 * ///
 * /// If targeted child is empty, place the particle there
 * ///
 * if (static_cast<cellPtrT>(curNodePtr)->child[targetOctant] == NULL)
 *  {
 *    curNodePtr->isEmpty = false;
 *
 *    newPartChild(targetOctant);
 *    goChild(targetOctant);
 *
 *    if (_newPartIdx >= noLocParts)
 *      {
 *        curNodePtr->isLocal = false;
 *      }
 *    else
 *      {
 *        curNodePtr->isLocal = true;
 *        /// save the particle's proxy
 *        partProxies[_newPartIdx] = curNodePtr;
 *      }
 *
 *    /// particle saves its position to node directly
 *    static_cast<partPtrT>(curNodePtr)->xPos = pos(_newPartIdx, X);
 *    static_cast<partPtrT>(curNodePtr)->yPos = pos(_newPartIdx, Y);
 *    static_cast<partPtrT>(curNodePtr)->zPos = pos(_newPartIdx, Z);
 *    static_cast<partPtrT>(curNodePtr)->mass = m(_newPartIdx);
 *
 *    /// ident saves the rowIndex of the particle
 *    curNodePtr->ident = _newPartIdx;
 *
 *    goUp();
 *  }
 *
 * ///
 * /// ... or if existing child is a node, then try to place the particle
 * /// as a child of this node
 * ///
 * else if (not
 *         static_cast<cellPtrT>(curNodePtr)->child[targetOctant]->isParticle)
 *  {
 *    curNodePtr->isEmpty = false;
 *    goChild(targetOctant);
 *    insertParticleRecursor(_newPartIdx);
 *    goUp();
 *  }
 *
 * ///
 * /// ... or if existing child is a particle (ghost/nonghost), then
 * /// replace it by a new node and try to place the existing two
 * /// particles as childs of this node
 * else if (static_cast<cellPtrT>(curNodePtr)->child[targetOctant]->isParticle)
 *  {
 *    /// goto child, save resident particle
 *    goChild(targetOctant);
 *    size_t residentIdx = static_cast<partPtrT>(curNodePtr)->ident;
 *    goUp();
 *
 *    /// replace particle by cell node
 *    delete static_cast<cellPtrT>(curNodePtr)->child[targetOctant];
 *    static_cast<cellPtrT>(curNodePtr)->child[targetOctant] = NULL;
 *    newCellChild(targetOctant);
 *
 *    /// and try to insert both particles again
 *    goChild(targetOctant);
 *    ///
 *    /// check whether we are too deep
 *    ///
 *    if (curNodePtr->depth > 128)
 *      throw PartsTooClose(curNodePtr->depth,
 *                          residentIdx,
 *                          _newPartIdx,
 *                          rootPtr);
 *
 *    insertParticleRecursor(residentIdx);
 *    insertParticleRecursor(_newPartIdx);
 *    goUp();
 *  }
 * else
 *  {
 *  }
 * }
 */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     // end of insertParticle() stuff

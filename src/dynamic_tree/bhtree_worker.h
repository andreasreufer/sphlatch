#ifndef BHTREE_WORKER_H
#define BHTREE_WORKER_H

/*
 *  bhtree_worker.h
 *
 *  Created by Andreas Reufer on 02.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_dynamic.h"
#include "bhtree_nodes.h"

#include "particle_manager.h"

namespace sphlatch {
class BHTreeWorker {
public:
   typedef BHTree&                treeRefT;
   typedef BHTree*                treePtrT;

   typedef ParticleManager        partManagerT;

   typedef BHTree::nodePtrT       nodePtrT;
   typedef BHTree::partPtrT       partPtrT;
   typedef BHTree::cellPtrT       cellPtrT;
   typedef BHTree::czllPtrT       czllPtrT;

   typedef BHTree::czllPtrListT   czllPtrListT;

   BHTreeWorker(const BHTreeWorker& _worker);
   BHTreeWorker(treePtrT _treePtr);
   ~BHTreeWorker();

protected:
   void goRoot();
   void goUp();
   void goNext();
   void goSkip();
   void goChild(const size_t _n);
   void goSibling(const size_t _n);

   bool pointInsideCell(const fType _x, const fType _y, const fType _z);
   size_t getOctant(const fType _x, const fType _y, const fType _z);

   nodePtrT curPtr, rootPtr;
   treePtrT treePtr;

public:
   partManagerT& PartManager;

   matrixRefType  pos, acc;
   valvectRefType eps, m, h, cost;

   void setTree(treePtrT _treePtr);
};

///
/// the copy constructor, mainly used where several workers
/// are spawned by OpenMP and need to be initialized from
/// the original worker
///
BHTreeWorker::BHTreeWorker(const BHTreeWorker& _worker) :
   PartManager(partManagerT::instance()),
   pos(PartManager.pos),
   acc(PartManager.acc),
   eps(PartManager.eps),
   m(PartManager.m),
   h(PartManager.h),
   cost(PartManager.cost)
{
   setTree(_worker.treePtr);
}

///
/// the standard constructor for a Tree worker
/// with the tree pointer as an argument
///
BHTreeWorker::BHTreeWorker(treePtrT _treePtr) :
   PartManager(partManagerT::instance()),
   pos(PartManager.pos),
   acc(PartManager.acc),
   eps(PartManager.eps),
   m(PartManager.m),
   h(PartManager.h),
   cost(PartManager.cost)
{
   setTree(_treePtr);
}

BHTreeWorker::~BHTreeWorker()
{ }

void BHTreeWorker::setTree(treePtrT _treePtr)
{
   treePtr = _treePtr;
   rootPtr = _treePtr->rootPtr;
   curPtr  = _treePtr->rootPtr;
}

void BHTreeWorker::goRoot()
{
   curPtr = rootPtr;
}

void BHTreeWorker::goUp()
{
   curPtr = curPtr->parent;
   assert(curPtr != NULL);
}

void BHTreeWorker::goNext()
{
   curPtr = curPtr->next;
   assert(curPtr != NULL);
}

void BHTreeWorker::goSkip()
{
   curPtr = curPtr->skip;
   assert(curPtr != NULL);
}

void BHTreeWorker::goChild(const size_t _n)
{
   curPtr = static_cast<cellPtrT>(curPtr)->child[_n];
   assert(curPtr != NULL);
}

void BHTreeWorker::goSibling(const size_t _n)
{
   curPtr = static_cast<czllPtrT>(curPtr)->sibling[_n];
   assert(curPtr != NULL);
}

bool BHTreeWorker::pointInsideCell(const fType _x,
                                   const fType _y,
                                   const fType _z)
{
   assert(curPtr != NULL);
   const fType hclSz = 0.5 * static_cast<cellPtrT>(curPtr)->clSz;

   // also try opposite logic for performance tests ...
   return(static_cast<cellPtrT>(curPtr)->xCen - hclSz < _x &&
          static_cast<cellPtrT>(curPtr)->xCen + hclSz > _x &&
          static_cast<cellPtrT>(curPtr)->yCen - hclSz < _y &&
          static_cast<cellPtrT>(curPtr)->yCen + hclSz > _y &&
          static_cast<cellPtrT>(curPtr)->zCen - hclSz < _z &&
          static_cast<cellPtrT>(curPtr)->zCen + hclSz > _z);
}

size_t BHTreeWorker::getOctant(const fType _x,
                               const fType _y,
                               const fType _z)
{
   size_t targetOctant = 0;

   assert(curPtr != NULL);
   targetOctant += _x < static_cast<cellPtrT>(curPtr)->xCen ? 0 : 1;
   targetOctant += _y < static_cast<cellPtrT>(curPtr)->xCen ? 0 : 2;
   targetOctant += _z < static_cast<cellPtrT>(curPtr)->xCen ? 0 : 4;
   return(targetOctant);
}
};
#endif

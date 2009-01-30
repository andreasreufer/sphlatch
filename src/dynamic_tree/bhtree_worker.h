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
#include "bhtree_node_cells.h"
#include "bhtree_node_particle.h"

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
   typedef BHTree::gcllPtrT       gcllPtrT;

   typedef BHTree::czllPtrListT   czllPtrListT;

   BHTreeWorker(const BHTreeWorker& _worker);
   BHTreeWorker(const treePtrT _treePtr);
   ~BHTreeWorker();

protected:
   void goRoot();
   void goUp();
   void goNext();
   void goSkip();
   void goChild(const size_t _n);
   void goNeighbour(const size_t _n);

   bool pointInsideCell(const fType _x, const fType _y, const fType _z);
   bool pointInsideCell(const fType _x, const fType _y, const fType _z,
                        nodePtrT _node);
   size_t getOctant(const fType _x, const fType _y, const fType _z);
   size_t getOctant(nodePtrT _nodePtr);

   const treePtrT treePtr;
   nodePtrT       curPtr, rootPtr;

public:
   partManagerT& PartManager;

   matrixRefType  pos, acc;
   valvectRefType eps, m, h, cost;
};

///
/// the copy constructor, mainly used where several workers
/// are spawned by OpenMP and need to be initialized from
/// the original worker
///
BHTreeWorker::BHTreeWorker(const BHTreeWorker& _worker) :
   treePtr(_worker.treePtr),
   rootPtr(_worker.rootPtr),
   PartManager(partManagerT::instance()),
   pos(PartManager.pos),
   acc(PartManager.acc),
   eps(PartManager.eps),
   m(PartManager.m),
   h(PartManager.h),
   cost(PartManager.cost)
{ }

///
/// the standard constructor for a Tree worker
/// with the tree pointer as an argument
///
BHTreeWorker::BHTreeWorker(const treePtrT _treePtr) :
   treePtr(_treePtr),
   rootPtr(_treePtr->rootPtr),
   PartManager(partManagerT::instance()),
   pos(PartManager.pos),
   acc(PartManager.acc),
   eps(PartManager.eps),
   m(PartManager.m),
   h(PartManager.h),
   cost(PartManager.cost)
{ }

BHTreeWorker::~BHTreeWorker()
{ }

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

void BHTreeWorker::goNeighbour(const size_t _n)
{
   curPtr = static_cast<czllPtrT>(curPtr)->neighbour[_n];
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

bool BHTreeWorker::pointInsideCell(const fType    _x,
                                   const fType    _y,
                                   const fType    _z,
                                   const nodePtrT _node)
{
   assert(curPtr != NULL);
   const fType hclSz = 0.5 * static_cast<gcllPtrT>(_node)->clSz;

   // also try opposite logic for performance tests ...
   return(static_cast<gcllPtrT>(_node)->xCen - hclSz < _x &&
          static_cast<gcllPtrT>(_node)->xCen + hclSz > _x &&
          static_cast<gcllPtrT>(_node)->yCen - hclSz < _y &&
          static_cast<gcllPtrT>(_node)->yCen + hclSz > _y &&
          static_cast<gcllPtrT>(_node)->zCen - hclSz < _z &&
          static_cast<gcllPtrT>(_node)->zCen + hclSz > _z);
}

size_t BHTreeWorker::getOctant(const fType _x,
                               const fType _y,
                               const fType _z)
{
   size_t targetOctant = 0;

   assert(curPtr != NULL);
   targetOctant += _x < static_cast<cellPtrT>(curPtr)->xCen ? 0 : 1;
   targetOctant += _y < static_cast<cellPtrT>(curPtr)->yCen ? 0 : 2;
   targetOctant += _z < static_cast<cellPtrT>(curPtr)->zCen ? 0 : 4;
   return(targetOctant);
}

size_t BHTreeWorker::getOctant(nodePtrT _nodePtr)
{
   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(curPtr)->child[i] == _nodePtr)
      {
         return(i);
      }
   }
   return(8);
}
};
#endif

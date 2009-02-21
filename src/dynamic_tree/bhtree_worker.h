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
   typedef BHTree&                    treeRefT;
   typedef BHTree*                    treePtrT;

   typedef ParticleManager            partManagerT;

   typedef BHTree::nodePtrT           nodePtrT;
   typedef BHTree::nodePtrCT          nodePtrCT;

   typedef BHTree::partT              partT;
   typedef BHTree::partPtrT           partPtrT;

   typedef BHTree::cellT              cellT;
   typedef BHTree::cellPtrT           cellPtrT;

   typedef BHTree::czllT              czllT;
   typedef BHTree::czllPtrT           czllPtrT;

   typedef BHTree::gcllPtrT           gcllPtrT;

   typedef BHTree::czllPtrListT       czllPtrListT;
   typedef BHTree::czllPtrListItrT    czllPtrListItrT;
   typedef BHTree::czllPtrListCItrT   czllPtrListCItrT;

   typedef BHTree::partAllocT         partAllocT;
   typedef BHTree::cellAllocT         cellAllocT;
   typedef BHTree::czllAllocT         czllAllocT;

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
   size_t getOctant(const fType _x, const fType _y, const fType _z,
                    const gcllPtrT _cellPtr);

   size_t getChildNo(nodePtrT _nodePtr);

   void czllToCell(nodePtrT _nodePtr);
   void cellToCZll(nodePtrT _nodePtr);
   void partToCell(nodePtrT _nodePtr, const size_t _oct);

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
   assert(curPtr != NULL);
   curPtr = curPtr->next;
}

void BHTreeWorker::goSkip()
{
   assert(curPtr != NULL);
   curPtr = static_cast<gcllPtrT>(curPtr)->skip;
}

void BHTreeWorker::goChild(const size_t _n)
{
   assert(curPtr != NULL);
   curPtr = static_cast<gcllPtrT>(curPtr)->child[_n];
}

void BHTreeWorker::goNeighbour(const size_t _n)
{
   assert(curPtr != NULL);
   curPtr = static_cast<czllPtrT>(curPtr)->neighbour[_n];
}

// \todo try const fType &_x
bool BHTreeWorker::pointInsideCell(const fType _x,
                                   const fType _y,
                                   const fType _z)
{
   assert(curPtr != NULL);
   const fType hclSz = 0.5 * static_cast<cellPtrT>(curPtr)->clSz;

   const fType xCen = static_cast<cellPtrT>(curPtr)->xCen;
   const fType yCen = static_cast<cellPtrT>(curPtr)->yCen;
   const fType zCen = static_cast<cellPtrT>(curPtr)->zCen;

   // \todo also try opposite logic for performance tests ...
   return(xCen - hclSz < _x && xCen + hclSz > _x &&
          yCen - hclSz < _y && yCen + hclSz > _y &&
          zCen - hclSz < _z && zCen + hclSz > _z);
}

bool BHTreeWorker::pointInsideCell(const fType    _x,
                                   const fType    _y,
                                   const fType    _z,
                                   const nodePtrT _node)
{
   assert(curPtr != NULL);
   const fType hclSz = 0.5 * static_cast<gcllPtrT>(_node)->clSz;

   const fType xCen = static_cast<cellPtrT>(_node)->xCen;
   const fType yCen = static_cast<cellPtrT>(_node)->yCen;
   const fType zCen = static_cast<cellPtrT>(_node)->zCen;

   // \todo also try opposite logic for performance tests ...
   return(xCen - hclSz < _x && xCen + hclSz > _x &&
          yCen - hclSz < _y && yCen + hclSz > _y &&
          zCen - hclSz < _z && zCen + hclSz > _z);
}

size_t BHTreeWorker::getOctant(const fType _x,
                               const fType _y,
                               const fType _z)
{
   size_t targetOctant = 0;

   assert(curPtr != NULL);
   targetOctant += _x < static_cast<gcllPtrT>(curPtr)->xCen ? 0 : 1;
   targetOctant += _y < static_cast<gcllPtrT>(curPtr)->yCen ? 0 : 2;
   targetOctant += _z < static_cast<gcllPtrT>(curPtr)->zCen ? 0 : 4;
   return(targetOctant);
}

size_t BHTreeWorker::getOctant(const fType    _x,
                               const fType    _y,
                               const fType    _z,
                               const gcllPtrT _cellPtr)
{
   size_t targetOctant = 0;

   assert(_cellPtr != NULL);
   targetOctant += _x < static_cast<gcllPtrT>(_cellPtr)->xCen ? 0 : 1;
   targetOctant += _y < static_cast<gcllPtrT>(_cellPtr)->yCen ? 0 : 2;
   targetOctant += _z < static_cast<gcllPtrT>(_cellPtr)->zCen ? 0 : 4;
   return(targetOctant);
}

size_t BHTreeWorker::getChildNo(nodePtrT _nodePtr)
{
   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(curPtr)->child[i] == _nodePtr)
         return(i);
   }
   return(8);
}

void BHTreeWorker::czllToCell(nodePtrT _nodePtr)
{
   const cellPtrT newCellPtr = treePtr->cellAllocator.pop();

   newCellPtr->initFromCZll(static_cast<czllT&>(*_nodePtr));
   treePtr->czllAllocator.push(static_cast<czllPtrT>(_nodePtr));
   _nodePtr = newCellPtr;
}

void BHTreeWorker::cellToCZll(nodePtrT _nodePtr)
{
   const czllPtrT newCZllPtr = treePtr->czllAllocator.pop();

   newCZllPtr->initFromCell(static_cast<cellT&>(*_nodePtr));
   treePtr->cellAllocator.push(static_cast<cellPtrT>(_nodePtr));
   _nodePtr = newCZllPtr;
}

void BHTreeWorker::partToCell(nodePtrT _cellPtr, const size_t _oct)
{
   assert(static_cast<gcllPtrT>(_cellPtr)->child[_oct] != NULL);
   assert(static_cast<gcllPtrT>(_cellPtr)->child[_oct]->isParticle);
   const partPtrT resPartPtr =
      static_cast<partPtrT>(static_cast<gcllPtrT>(_cellPtr)->child[_oct]);

   const cellPtrT newCellPtr = treePtr->cellAllocator.pop();

   newCellPtr->clear();
   newCellPtr->ident = treePtr->noCells;
   treePtr->noCells++;

   // wire new cell and its parent
   newCellPtr->parent = _cellPtr;
   static_cast<gcllPtrT>(_cellPtr)->child[_oct] = newCellPtr;

   // set cell position
   newCellPtr->inheritCellPos(_oct);

   // add cost of the particle to be added
   static_cast<gcllPtrT>(newCellPtr)->cost    = resPartPtr->cost;
   static_cast<gcllPtrT>(newCellPtr)->noParts = 1;

   // wire particle to new node
   const fType  posX   = resPartPtr->xPos;
   const fType  posY   = resPartPtr->yPos;
   const fType  posZ   = resPartPtr->zPos;
   const size_t newOct = getOctant(posX, posY, posZ, newCellPtr);

   newCellPtr->child[newOct] = resPartPtr;
   resPartPtr->parent        = newCellPtr;
}

///
/// specialization for a read-only BHTree worker (eg. gravity walk)
///
class BHTreeWorkerRO : public BHTreeWorker {
public:
   BHTreeWorkerRO(treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeWorkerRO(const BHTreeWorkerRO& _workerRO) : BHTreeWorker(_workerRO) { }

protected:
   ///
   /// overwrite curPtr by a faster read-only version
   ///
   nodePtrCT curPtr, rootPtr;
};
};
#endif

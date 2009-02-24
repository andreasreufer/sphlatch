#ifndef BHTREE_WORKER_H
#define BHTREE_WORKER_H

/*
 *  bhtree_worker.h
 *
 *  Created by Andreas Reufer on 02.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#ifdef SPHLATCH_OPENMP
 #include <omp.h>
#endif

#include "bhtree_dynamic.h"
#include "bhtree_node_cells.h"
#include "bhtree_node_particle.h"

namespace sphlatch {
class BHTreeWorker {
public:
   typedef BHTree&                    treeRefT;
   typedef BHTree*                    treePtrT;

   typedef BHTree::nodePtrT           nodePtrT;
   typedef BHTree::nodePtrCT          nodePtrCT;

   typedef BHTree::partT              partT;
   typedef BHTree::partPtrT           partPtrT;

   typedef BHTree::pnodT              pnodT;
   typedef BHTree::pnodPtrT           pnodPtrT;

   typedef BHTree::cellT              cellT;
   typedef BHTree::cellPtrT           cellPtrT;

   typedef BHTree::czllT              czllT;
   typedef BHTree::czllPtrT           czllPtrT;

   typedef BHTree::gcllPtrT           gcllPtrT;

   typedef BHTree::czllPtrListT       czllPtrListT;
   typedef BHTree::czllPtrListItrT    czllPtrListItrT;
   typedef BHTree::czllPtrListCItrT   czllPtrListCItrT;

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

   bool pointInsideCell(const vect3dT& _pos);
   bool pointInsideCell(const vect3dT& _pos, const nodePtrT _nodePtr);

   size_t getOctant(const vect3dT& _pos);
   size_t getOctant(const vect3dT& _pos, const nodePtrT _nodePtr);

   size_t getChildNo(nodePtrT _nodePtr);

   void czllToCell(nodePtrT _nodePtr);
   void cellToCZll(nodePtrT _nodePtr);
   void partToCell(nodePtrT _nodePtr, const size_t _oct);

   const size_t   noThreads, myThread;
   const treePtrT treePtr;

   nodePtrT curPtr, rootPtr;
};

///
/// the copy constructor, mainly used where several workers
/// are spawned by OpenMP and need to be initialized from
/// the original worker
///
BHTreeWorker::BHTreeWorker(const BHTreeWorker& _worker) :
#ifdef SPHLATCH_OPENMP
   noThreads(omp_get_num_threads()),
   myThread(omp_get_thread_num()),
#else
   noThreads(1),
   myThread(0),
#endif
   treePtr(_worker.treePtr),
   rootPtr(_worker.rootPtr) //,
{ }

///
/// the standard constructor for a Tree worker
/// with the tree pointer as an argument
///
BHTreeWorker::BHTreeWorker(const treePtrT _treePtr) :
#ifdef SPHLATCH_OPENMP
   noThreads(omp_get_num_threads()),
   myThread(omp_get_thread_num()),
#else
   noThreads(1),
   myThread(0),
#endif
   treePtr(_treePtr),
   rootPtr(_treePtr->rootPtr)
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

bool BHTreeWorker::pointInsideCell(const vect3dT& _pos)
{
   assert(curPtr != NULL);
   const fType hclSz = 0.5 * static_cast<cellPtrT>(curPtr)->clSz;

   return(all(static_cast<gcllPtrT>(curPtr)->cen - _pos < hclSz) &&
          all(static_cast<gcllPtrT>(curPtr)->cen - _pos > -hclSz));
}

bool BHTreeWorker::pointInsideCell(const vect3dT& _pos, nodePtrT _nodePtr)
{
   assert(_nodePtr != NULL);
   const fType hclSz = 0.5 * static_cast<cellPtrT>(_nodePtr)->clSz;

   return(all(static_cast<gcllPtrT>(_nodePtr)->cen - _pos < hclSz) &&
          all(static_cast<gcllPtrT>(_nodePtr)->cen - _pos > -hclSz));
}

size_t BHTreeWorker::getOctant(const vect3dT& _pos)
{
   size_t targetOctant = 0;

   assert(_cellPtr != NULL);
   targetOctant += _pos[0] < static_cast<gcllPtrT>(curPtr)->cen[0] ? 0 : 1;
   targetOctant += _pos[1] < static_cast<gcllPtrT>(curPtr)->cen[1] ? 0 : 2;
   targetOctant += _pos[2] < static_cast<gcllPtrT>(curPtr)->cen[2] ? 0 : 4;
   return(targetOctant);
}

size_t BHTreeWorker::getOctant(const vect3dT& _pos, const nodePtrT _nodePtr)
{
   size_t targetOctant = 0;

   assert(_cellPtr != NULL);
   targetOctant += _pos[0] < static_cast<gcllPtrT>(_nodePtr)->cen[0] ? 0 : 1;
   targetOctant += _pos[1] < static_cast<gcllPtrT>(_nodePtr)->cen[1] ? 0 : 2;
   targetOctant += _pos[2] < static_cast<gcllPtrT>(_nodePtr)->cen[2] ? 0 : 4;
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
   const cellPtrT newCellPtr = new cellT;

   newCellPtr->initFromCZll(static_cast<czllT&>(*_nodePtr));
   delete static_cast<czllPtrT>(_nodePtr);
   _nodePtr = newCellPtr;
}

void BHTreeWorker::cellToCZll(nodePtrT _nodePtr)
{
   const czllPtrT newCZllPtr = new czllT;

   newCZllPtr->initFromCell(static_cast<cellT&>(*_nodePtr));
   delete static_cast<cellPtrT>(_nodePtr);
   _nodePtr = newCZllPtr;
}

void BHTreeWorker::partToCell(nodePtrT _cellPtr, const size_t _oct)
{
   assert(static_cast<gcllPtrT>(_cellPtr)->child[_oct] != NULL);
   assert(static_cast<gcllPtrT>(_cellPtr)->child[_oct]->isParticle);
   const pnodPtrT resPartPtr =
      static_cast<pnodPtrT>(static_cast<gcllPtrT>(_cellPtr)->child[_oct]);

   const cellPtrT newCellPtr = new cellT;

   newCellPtr->clear();
   newCellPtr->ident = treePtr->noCells;
   treePtr->noCells++;

   // wire new cell and its parent
   newCellPtr->parent = _cellPtr;
   static_cast<gcllPtrT>(_cellPtr)->child[_oct] = newCellPtr;

   // set cell position
   newCellPtr->inheritCellPos(_oct);

   // add cost of the particle to be added
   static_cast<gcllPtrT>(newCellPtr)->cost    = resPartPtr->partPtr->cost;
   static_cast<gcllPtrT>(newCellPtr)->noParts = 1;

   // wire particle to new node
   const size_t newOct = getOctant(resPartPtr->pos, newCellPtr);

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

#ifndef BHTREE_WORKER_CPP
#define BHTREE_WORKER_CPP

/*
 *  bhtree_worker.cpp
 *
 *  Created by Andreas Reufer on 07.10.09
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#ifdef SPHLATCH_OPENMP
 #include <omp.h>
#endif

#include "bhtree_worker.h"
#include "bhtree.cpp"

namespace sphlatch {
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
   rootPtr(_worker.rootPtr)
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
   const fType hclSz = 0.5 * static_cast<gcllPtrT>(curPtr)->clSz;

   return(all(static_cast<gcllPtrT>(curPtr)->cen - _pos < hclSz) &&
          all(static_cast<gcllPtrT>(curPtr)->cen - _pos > -hclSz));
}

bool BHTreeWorker::pointInsideCell(const vect3dT& _pos, nodePtrT _nodePtr)
{
   assert(_nodePtr != NULL);
   const fType hclSz = 0.5 * static_cast<gcllPtrT>(_nodePtr)->clSz;

   return(all(static_cast<gcllPtrT>(_nodePtr)->cen - _pos < hclSz) &&
          all(static_cast<gcllPtrT>(_nodePtr)->cen - _pos > -hclSz));
}

bool BHTreeWorker::sphereTotInCell(const vect3dT& _pos, const fType& _r)
{
   assert(curPtr != NULL);
   const fType RCellMRSph = 0.5 * static_cast<gcllPtrT>(curPtr)->clSz - _r;
   return(all(static_cast<gcllPtrT>(curPtr)->cen - _pos < RCellMRSph) &&
          all(static_cast<gcllPtrT>(curPtr)->cen - _pos > -RCellMRSph));
}

bool BHTreeWorker::sphereTotOutCell(const vect3dT& _pos, const fType& _r)
{
   assert(curPtr != NULL);
   const fType RCellPRSph = 0.5 * static_cast<gcllPtrT>(curPtr)->clSz + _r;

   return(any(static_cast<gcllPtrT>(curPtr)->cen - _pos > RCellPRSph) ||
          any(static_cast<gcllPtrT>(curPtr)->cen - _pos < -RCellPRSph));
}

///
/// find maximal mass enclosing radius:
/// - load current particle i data
/// - go to current particle i node
/// - go up until the cell j contains enough mass m
/// - the sphere around the particle with r_s = | r_j - r_i + 0.5*r_c |
///   contains at least enough mass m
///   r_c is the distance between the cell center and the cell
///   corner with the biggest distance to the particle
///
/// this algorithm gives a maximal smoothing length for each particle:
/// when each particle has mass m = 1, then the masses of the cells
/// give directly the number of particles contained in them
///
fType BHTreeWorker::maxMassEncloseRad(const pnodPtrT _partPtr, const fType _m)
{
   assert(_partPtr != NULL);
   curPtr = _partPtr;

   const vect3dT ppos = static_cast<pnodPtrT>(curPtr)->pos;

   while (curPtr->parent != NULL)
   {
      goUp();
      if (static_cast<mcllPtrT>(curPtr)->m > _m)
         break;
   }

   const fType hcSize = 0.5 * static_cast<gcllPtrT>(curPtr)->clSz;

   // search for the farest corner from the particle
   vect3dT corner = static_cast<gcllPtrT>(curPtr)->cen;
   for (size_t i = 0; i < 3; i++)
      corner[i] = ppos[i] < corner[i] ? corner[i] + hcSize : corner[i] - hcSize;

   return(sqrt(dot(ppos - corner, ppos - corner)));
}

///
/// get the cell octant for a position, no checks are performed whether
/// the position is actually in the cell
///
size_t BHTreeWorker::getOctant(const vect3dT& _pos)
{
   size_t targetOctant = 0;

   assert(curPtr != NULL);
   targetOctant += _pos[0] < static_cast<gcllPtrT>(curPtr)->cen[0] ? 0 : 1;
   targetOctant += _pos[1] < static_cast<gcllPtrT>(curPtr)->cen[1] ? 0 : 2;
   targetOctant += _pos[2] < static_cast<gcllPtrT>(curPtr)->cen[2] ? 0 : 4;
   return(targetOctant);
}

///
/// get the cell octant for a position, no checks are performed whether
/// the position is actually in the cell
///
size_t BHTreeWorker::getOctant(const vect3dT& _pos, const nodePtrT _nodePtr)
{
   size_t targetOctant = 0;

   assert(_nodePtr != NULL);
   targetOctant += _pos[0] < static_cast<gcllPtrT>(_nodePtr)->cen[0] ? 0 : 1;
   targetOctant += _pos[1] < static_cast<gcllPtrT>(_nodePtr)->cen[1] ? 0 : 2;
   targetOctant += _pos[2] < static_cast<gcllPtrT>(_nodePtr)->cen[2] ? 0 : 4;
   return(targetOctant);
}

///
/// get the child index for a node ptr
///
size_t BHTreeWorker::getChildNo(nodePtrT _nodePtr)
{
   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(curPtr)->child[i] == _nodePtr)
         return(i);
   }
   return(8);
}

size_t BHTreeWorker::getChildNo(nodePtrT _nodePtr, nodePtrT _parPtr)
{
   for (size_t i = 0; i < 8; i++)
   {
      if (static_cast<gcllPtrT>(_parPtr)->child[i] == _nodePtr)
         return(i);
   }
   return(8);
}

///
/// transforms a CZ cell into a quadrupole
///
qcllPtrT BHTreeWorker::czllToCell(nodePtrT _nodePtr)
{
   const qcllPtrT newCellPtr = new qcllT;

   newCellPtr->initFromCZll(static_cast<czllT&>(*_nodePtr));
   delete static_cast<czllPtrT>(_nodePtr);
   return(newCellPtr);
}

///
/// transforms a quadrupole cell into a CZ cell
///
czllPtrT BHTreeWorker::cellToCZll(nodePtrT _nodePtr)
{
   const czllPtrT newCZllPtr = new czllT;

   newCZllPtr->initFromCell(static_cast<qcllT&>(*_nodePtr));
   delete static_cast<qcllPtrT>(_nodePtr);
   return(newCZllPtr);
}

qcllPtrT BHTreeWorker::partToCell(nodePtrT _cellPtr, const size_t _oct)
{
   assert(static_cast<gcllPtrT>(_cellPtr)->child[_oct] != NULL);
   assert(static_cast<gcllPtrT>(_cellPtr)->child[_oct]->isParticle);
   const pnodPtrT resPartPtr =
      static_cast<pnodPtrT>(static_cast<gcllPtrT>(_cellPtr)->child[_oct]);

   const qcllPtrT newCellPtr = new qcllT;

   newCellPtr->clear();
   newCellPtr->ident = treePtr->noCells;
   treePtr->noCells++;

   // wire new cell and its parent
   newCellPtr->parent = _cellPtr;
   static_cast<gcllPtrT>(_cellPtr)->child[_oct] = newCellPtr;

   // set cell position
   newCellPtr->inheritCellPos(_oct);

   // wire particle to new node
   const size_t newOct = getOctant(resPartPtr->pos, newCellPtr);

   newCellPtr->child[newOct] = resPartPtr;
   resPartPtr->parent        = newCellPtr;
   resPartPtr->depth         = newCellPtr->depth + 1;

   return(newCellPtr);
}
};
#endif

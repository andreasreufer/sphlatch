#ifndef BHTREE_WORKER_SPHSUM_CPP
#define BHTREE_WORKER_SPHSUM_CPP

/*
 *  bhtree_worker_sphsum.cpp
 *
 *  Created by Andreas Reufer on 14.10.09.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"
#include "bhtree_particle.h"
#include "timer.cpp"

#include "bhtree_worker_neighfunc.cpp"

namespace sphlatch {
template<typename _sumT, typename _partT>
class SPHsumWorker : public BHTreeWorker {
public:
   SPHsumWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   SPHsumWorker(const SPHsumWorker& _SPHwork) : BHTreeWorker(_SPHwork) { }
   ~SPHsumWorker() { }

   void operator()(const czllPtrT _czll);

   typedef sphlatch::Timer   timerT;

private:
   _sumT  Sum;
   timerT Timer;
   void sumNeighbours(const pnodPtrT _part);
};

template<typename _sumT, typename _partT>
void SPHsumWorker<_sumT, _partT>::operator()(const czllPtrT _czll)
{
   nodePtrT       curPart  = _czll->chldFrst;
   const nodePtrT stopChld = _czll->chldLast->next;

   // an empty CZ cell may have an chldFrst pointing to NULL
   if (curPart == NULL)
      return;


   Timer.start();
   while (curPart != stopChld)
   {
      if (curPart->isParticle)
         sumNeighbours(static_cast<pnodPtrT>(curPart));
      curPart = curPart->next;
   }
   const double compTime = Timer.getRoundTime();
   _czll->compTime += static_cast<fType>(compTime);
}

template<typename _sumT, typename _partT>
void SPHsumWorker<_sumT, _partT>::sumNeighbours(const pnodPtrT _part)
{
   // go to the particle and load its data
   curPtr = _part;
   _partT* const ipartPtr = static_cast<_partT*>(_part->partPtr);

   Sum.preSum(ipartPtr);

   const vect3dT ppos = _part->pos;

   //FIXME: this factor should be set more generically
   const fType hi    = static_cast<_partT*>(_part->partPtr)->h;
   const fType srad  = 2. * hi;
   const fType srad2 = srad * srad;

#ifdef SPHLATCH_NONEIGH
   size_t non = 0;
#endif

   // go to particles parent cell
   goUp();

   // go up, until the search sphere is completely in the current cell
   while (not sphereTotInCell(ppos, srad) && curPtr->parent != NULL)
      goUp();


   // now start to search the subtree for potential neighbours
   const nodePtrT lastNode = static_cast<gcllPtrT>(curPtr)->skip;
   while (curPtr != lastNode)
   {
      if (curPtr->isParticle)
      {
         const vect3dT rvec = ppos - static_cast<pnodPtrT>(curPtr)->pos;
         const fType   rr   = dot(rvec, rvec);

         if (rr < srad2)
         {
            Sum(ipartPtr,
                static_cast<_partT*>(static_cast<pnodPtrT>(curPtr)->partPtr),
                rvec, rr, srad);
#ifdef SPHLATCH_NONEIGH
            non++;
#endif
         }

         goNext();
      }
      else
      {
         // if search sphere completely outside of the current cell, skip it
         if (sphereTotOutCell(ppos, srad))
         {
            if (static_cast<gcllPtrT>(curPtr)->skip == NULL)
               break;
            else
               goSkip();
         }
         else
            goNext();
      }
   }
#ifdef SPHLATCH_NONEIGH
   static_cast<_partT*>(_part->partPtr)->noneigh = non;
#endif

   Sum.postSum(ipartPtr);
}

template<typename _sumT, typename _partT>
class SPHsum2Worker : public NeighWorker<_sumT, _partT> {
public:
   SPHsum2Worker(const BHTreeWorker::treePtrT _treePtr) :
      NeighWorker<_sumT, _partT>(_treePtr) { }
   SPHsum2Worker(const SPHsum2Worker& _SPHwork) :
      NeighWorker<_sumT, _partT>(_SPHwork) { }
   ~SPHsum2Worker() { }

   void operator()(const czllPtrT _czll);

   typedef sphlatch::Timer   timerT;
private:

   _sumT  Sum;
   timerT Timer;
};


template<typename _sumT, typename _partT>
void SPHsum2Worker<_sumT, _partT>::operator()(const czllPtrT _czll)
{
   nodePtrT       curPart  = _czll->chldFrst;
   const nodePtrT stopChld = _czll->chldLast->next;

   // an empty CZ cell may have an chldFrst pointing to NULL
   if (curPart == NULL)
      return;

   Timer.start();
   while (curPart != stopChld)
   {
      if (curPart->isParticle)
      {
         _partT* const partPtr = static_cast<_partT*>(
           static_cast<pnodPtrT>(curPart)->partPtr);
         const fType hi = partPtr->h;
         const fType srad = 2. * hi;
         
         Sum.preSum(partPtr);
         NeighWorker<_sumT,
                     _partT>::neighExecFunc(static_cast<pnodPtrT>(curPart),
                                            srad);
         Sum.postSum(partPtr);
      }
      curPart = curPart->next;
   }
   const double compTime = Timer.getRoundTime();
   _czll->compTime += static_cast<fType>(compTime);
}
};


#endif

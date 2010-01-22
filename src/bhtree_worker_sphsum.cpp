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
class SPHsumWorker : public NeighWorker<_sumT, _partT> {
public:
   SPHsumWorker(const BHTreeWorker::treePtrT _treePtr) :
      NeighWorker<_sumT, _partT>(_treePtr)
   { }

   SPHsumWorker(const SPHsumWorker& _SPHwork) :
      NeighWorker<_sumT, _partT>(_SPHwork)
   { }

   ~SPHsumWorker()
   { }

   void operator()(const czllPtrT _czll);
   void operator()(_partT* const _partPtr);

   typedef sphlatch::Timer   timerT;

private:
   timerT Timer;
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
      {
         _partT* const partPtr = static_cast<_partT*>(
            static_cast<pnodPtrT>(curPart)->partPtr);
         const fType hi   = partPtr->h;
         const fType srad = 2. * hi;

         NeighWorker<_sumT, _partT>::Func.preSum(partPtr);
         NeighWorker<_sumT,
                     _partT>::neighExecFunc(static_cast<pnodPtrT>(curPart),
                                            srad);
         NeighWorker<_sumT, _partT>::Func.postSum(partPtr);
      }
      curPart = curPart->next;
   }
   const double compTime = Timer.getRoundTime();
   _czll->compTime += static_cast<fType>(compTime);
}

template<typename _sumT, typename _partT>
void SPHsumWorker<_sumT, _partT>::operator()(_partT* const _partPtr)
{
   pnodPtrT const pnodPtr = _partPtr->treeNode;

   const fType hi   = _partPtr->h;
   const fType srad = 2. * hi;

   NeighWorker<_sumT, _partT>::Func.preSum(_partPtr);
   NeighWorker<_sumT, _partT>::neighExecFunc(static_cast<pnodPtrT>(pnodPtr),
                                             srad);
   NeighWorker<_sumT, _partT>::Func.postSum(_partPtr);
}
};
#endif

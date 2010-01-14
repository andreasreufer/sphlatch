#ifndef BHTREE_WORKER_COST_CPP
#define BHTREE_WORKER_COST_CPP

/*
 *  bhtree_worker_cost.cpp
 *
 *  Created by Andreas Reufer on 14.01.10
 *  Copyright 2010 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"
#include "bhtree.h"

namespace sphlatch {
template<typename _partT>
class CostWorker : public BHTreeWorker {
public:
   CostWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   CostWorker(const CostWorker& _cw) : BHTreeWorker(_cw) { }
   ~CostWorker() { }

   void operator()(const czllPtrT _czll);

private:
};

template<typename _partT>
void CostWorker<_partT>::operator()(const czllPtrT _czll)
{
   nodePtrT       curPart  = _czll->chldFrst;
   const nodePtrT stopChld = _czll->chldLast->next;

   const fType totTime = _czll->compTime;
   const fType absTime = totTime / static_cast<fType>(_czll->noParts);

   // an empty CZ cell may have an chldFrst pointing to NULL
   if (curPart == NULL)
      return;

   while (curPart != stopChld)
   {
      if (curPart->isParticle)
        static_cast<pnodPtrT>(curPart)->partPtr->cost = absTime;
      curPart = curPart->next;
   }
}

};

#endif

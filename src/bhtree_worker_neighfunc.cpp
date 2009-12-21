#ifndef BHTREE_WORKER_NEIGHFUNC_CPP
#define BHTREE_WORKER_NEIGHFUNC_CPP

/*
 *  bhtree_worker_neighfunc.cpp
 *
 *  Created by Andreas Reufer on 14.10.09.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"
#include "bhtree_particle.h"
#include "timer.cpp"

namespace sphlatch {
template<typename _funcT, typename _partT>
class NeighWorker : public BHTreeWorker {
public:
   NeighWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   NeighWorker(const NeighWorker& _Nworker) : BHTreeWorker(_Nworker) { }
   ~NeighWorker() { }

   void neighExecFunc(const pnodPtrT _part, const fType _srad);
   void neighExecFunc(const czllPtrT _czll, const fType _srad);

private:
   _funcT Func;
};

template<typename _funcT, typename _partT>
void NeighWorker<_funcT, _partT>::neighExecFunc(const pnodPtrT _pnod,
                                                const fType    _srad)
{
   // go to the particle and load its data
   curPtr = _pnod;
   _partT* const ipartPtr = static_cast<_partT*>(_pnod->partPtr);
   const vect3dT ppos     = _pnod->pos;
   const fType srad2 = _srad*_srad;

   // go to particles parent cell
   goUp();

   // go up, until the search sphere is completely in the current cell
   while (not sphereTotInCell(ppos, _srad) && curPtr->parent != NULL)
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
            Func(ipartPtr,
                 static_cast<_partT*>(static_cast<pnodPtrT>(curPtr)->partPtr),
                 rvec, rr, _srad);
         }

         goNext();
      }
      else
      {
         // if search sphere completely outside of the current cell, skip it
         if (sphereTotOutCell(ppos, _srad))
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
}

template<typename _funcT, typename _partT>
void NeighWorker<_funcT, _partT>::neighExecFunc(const czllPtrT _czll,
                                                const fType    _srad)
{
   nodePtrT       curPart  = _czll->chldFrst;
   const nodePtrT stopChld = _czll->chldLast->next;


   // an empty CZ cell may have an chldFrst pointing to NULL
   if (curPart == NULL)
      return;

   while (curPart != stopChld)
   {
     if (curPart->isParticle)
       neighExecFunc(static_cast<pnodPtrT>(curPart), _srad);
     curPart = curPart->next;
   }
}


};


#endif

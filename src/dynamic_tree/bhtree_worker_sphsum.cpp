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

namespace sphlatch {
template<typename _sumT, typename _partT>
class SPHsumWorker : public BHTreeWorker {
public:
   SPHsumWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   SPHsumWorker(const SPHsumWorker& _SPHwork) : BHTreeWorker(_SPHwork) { }
   ~SPHsumWorker() { }

   _sumT Sum;
   void operator()(const czllPtrT _czll);

private:
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

   //double time = omp_get_wtime();

   while (curPart != stopChld)
   {
      if (curPart->isParticle)
         sumNeighbours(static_cast<pnodPtrT>(curPart));
      curPart = curPart->next;
   }

   //std::cout << this << ":" << _czll << ":" << myThread << "   " << omp_get_wtime() - time << "\n";
}

template<typename _sumT, typename _partT>
void SPHsumWorker<_sumT, _partT>::sumNeighbours(const pnodPtrT _part)
{
   // go to the particle and load its data
   curPtr = _part;
   _partT* const ipartPtr = static_cast<_partT*>(_part->partPtr);
   const vect3dT ppos     = _part->pos;

   //FIXME: this factor should be set more generically
   const fType hi = static_cast<_partT*>(_part->partPtr)->h;
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

   Sum.preSum(ipartPtr);

   // now start to search the subtree for potential neighbours
   const nodePtrT lastNode = static_cast<gcllPtrT>(curPtr)->skip;
   while (curPtr != lastNode)
   {
      if (curPtr->isParticle)
      {
         const vect3dT rvec = ppos - static_cast<pnodPtrT>(curPtr)->pos;
         const fType rr = dot( rvec, rvec );

         if (rr < srad2)
         {
            Sum(ipartPtr,
                static_cast<_partT*>(static_cast<pnodPtrT>(curPtr)->partPtr),
                rvec, rr, hi);
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
};


#endif

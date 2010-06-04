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
#include <list>

namespace sphlatch {
template<typename _funcT, typename _partT>
class NeighWorker : public BHTreeWorker {
public:
   NeighWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   NeighWorker(const NeighWorker& _Nworker) : BHTreeWorker(_Nworker) { }
   ~NeighWorker() { }

   void neighExecFunc(const pnodPtrT _part, const fType _srad);
   void neighExecFunc(const czllPtrT _czll, const fType _srad);

protected:
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
   const fType   srad2    = _srad * _srad;

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

template<typename _partT>
class NeighListFunc
{
public:
   void operator()(_partT* const _i,
                   _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      neighList.push_back(_j);
   }

   std::list<_partT*> neighList;
};



template<typename _partT>
class NeighFindWorker : public NeighWorker<NeighListFunc<_partT> , _partT> {
public:
   typedef NeighWorker<NeighListFunc<_partT> , _partT>   parentT;

   NeighFindWorker(const BHTreeWorker::treePtrT _treePtr) :
      NeighWorker<NeighListFunc<_partT> , _partT>(_treePtr) { }
   NeighFindWorker(const NeighFindWorker& _NFwork) :
      NeighWorker<NeighListFunc<_partT> , _partT>(_NFwork) { }
   ~NeighFindWorker() { }

   std::list<_partT*> operator()(const _partT* _part, const fType _srad)
   {
      parentT::Func.neighList.clear();
      parentT::neighExecFunc(_part->treeNode, _srad);
      return(parentT::Func.neighList);
   }

   std::list<_partT*> operator()(const pnodPtrT _pnod, const fType _srad)
   {
      parentT::Func.neighList.clear();
      parentT::neighExecFunc(_pnod, _srad);
      return(parentT::Func.neighList);
   }
};

//
// small class which makes a list of neighbours and also store
// the squared distance to those neighbours
//
template<typename _partT>
class NeighDistListFunc
{
public:
   // small helper class which keeps a particle pointer and the squared dist.
   class PPtrD
   {
public:
      PPtrD(_partT* _ptr, const fType _rr) : ptr(_ptr), rr(_rr) { }
      ~PPtrD() { }
      _partT      * ptr;
      const fType rr;
   };

   // small helper class to sort a list of neighbours according to distance
   class PPtrDsorter
   {
public:
      bool operator()(PPtrD& _ppd1, PPtrD& _ppd2)
      {
         return(_ppd1.rr < _ppd2.rr);
      }
   };

   void operator()(_partT* const _i,
                   _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      neighList.push_back(PPtrD(_j, _rr));
   }

   std::list<PPtrD> neighList;
};


template<typename _partT>
class SmoLenFindWorker :
   public NeighWorker<NeighDistListFunc<_partT> , _partT> {
public:
   typedef NeighWorker<NeighDistListFunc<_partT> , _partT>   parentT;
   typedef typename NeighDistListFunc<_partT>::PPtrD         pptrDT;
   typedef typename NeighDistListFunc<_partT>::PPtrDsorter   pptrDsortT;
   typedef typename std::list<pptrDT>                        pptrDLT;

   SmoLenFindWorker(const BHTreeWorker::treePtrT _treePtr, const fType _mult) :
      NeighWorker<NeighDistListFunc<_partT> , _partT>(_treePtr), mult(_mult) { }
   SmoLenFindWorker(const SmoLenFindWorker& _SLFwork) :
      NeighWorker<NeighDistListFunc<_partT> , _partT>(_SLFwork)
      , mult(_SLFwork.mult) { }
   ~SmoLenFindWorker() { }

   void operator()(const czllPtrT _czll)
   {
      nodePtrT       curPart  = _czll->chldFrst;
      const nodePtrT stopChld = _czll->chldLast->next;

      if (curPart == NULL)
         return;

      while (curPart != stopChld)
      {
         if (curPart->isParticle)
         {
            const pnodPtrT pnod    = static_cast<pnodPtrT>(curPart);
            _partT* const  partPtr = static_cast<_partT*>(pnod->partPtr);

            const size_t noneigh = partPtr->noneigh;
            const fType  smass   = static_cast<fType>(noneigh);
            const fType  srad    = BHTreeWorker::maxMassEncloseRad(pnod, smass);

            assert(noneigh > 0);

            // clear the neighbour list and start the search
            parentT::Func.neighList.clear();
            parentT::neighExecFunc(pnod, srad);

            // sort the list according to the distance
            pptrDLT&   neighList(parentT::Func.neighList);
            pptrDsortT neighSorter;
            neighList.sort(neighSorter);

            if (neighList.size() <= noneigh)
               partPtr->h = sqrt((*(neighList.end())).rr) / mult;
            else
            {
               typename pptrDLT::const_iterator nitr = neighList.begin();
               for (size_t i = 0; i < noneigh; i++)
                  nitr++;

               partPtr->h = sqrt((*nitr).rr) / mult;
            }
         }
         curPart = curPart->next;
      }
   }

protected:
   const fType mult;
};


template<typename _partT>
class NeighFindSortedWorker :
   public NeighWorker<NeighDistListFunc<_partT> , _partT> {
public:
   typedef NeighWorker<NeighDistListFunc<_partT> , _partT>   parentT;
   typedef typename NeighDistListFunc<_partT>::PPtrD         pptrDT;
   typedef typename NeighDistListFunc<_partT>::PPtrDsorter   pptrDsortT;
   typedef typename std::list<pptrDT>                        pptrDLT;

   NeighFindSortedWorker(const BHTreeWorker::treePtrT _treePtr) :
      NeighWorker<NeighDistListFunc<_partT> , _partT>(_treePtr) { }
   NeighFindSortedWorker(const NeighFindSortedWorker& _SLFwork) :
      NeighWorker<NeighDistListFunc<_partT> , _partT>(_SLFwork) { }
   ~NeighFindSortedWorker() { }

   pptrDLT operator()(const _partT* _part, const fType _srad)
   {
      parentT::Func.neighList.clear();
      parentT::neighExecFunc(_part->treeNode, _srad);


      // sort the list according to the distance
      pptrDLT&   neighList(parentT::Func.neighList);
      pptrDsortT neighSorter;
      neighList.sort(neighSorter);

      return(parentT::Func.neighList);
   }
};
};

#endif

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

template<typename _partT>
class NeighDistListFunc
{
public:
   class PartPtrDist 
   {
     public:
       PartPtrDist(_partT* _ptr, const fType _rr) : ptr(_ptr), rr(_rr) {};
       ~PartPtrDist() {};
       _partT* ptr;
       const fType rr;
   };

   void operator()(_partT* const _i,
                   _partT* const _j,
                   const vect3dT& _rvec,
                   const fType _rr,
                   const fType _srad)
   {
      neighList.push_back( PartPtrDist(_j, _rr) );
   }

   std::list<PartPtrDist> neighList;
};

template<typename _partT>

class SmoLenFindWorker : public NeighWorker<NeighDistListFunc<_partT> , _partT> {
public:
   typedef NeighWorker<NeighDistListFunc<_partT> , _partT>   parentT;
   //typedef NeighDistListFunc<_partT>::PartPtrDist PartPtrDistT;
   
   //std::list<PartPtrDistT> neighList;

   SmoLenFindWorker(const BHTreeWorker::treePtrT _treePtr) :
      NeighWorker<NeighDistListFunc<_partT> , _partT>(_treePtr) { }
   SmoLenFindWorker(const SmoLenFindWorker& _NFwork) :
      NeighWorker<NeighDistListFunc<_partT> , _partT>(_NFwork) { }
   ~SmoLenFindWorker() { }

/*class SmoLenFindWorker : public NeighFindWorker<_partT> {
public:
   typedef NeighFindWorker<_partT>   parentT;

   SmoLenFindWorker(const BHTreeWorker::treePtrT _treePtr) :
      NeighFindWorker<_partT>(_treePtr) { }
   SmoLenFindWorker(const SmoLenFindWorker& _SLwork) :
      NeighFindWorker<_partT>(_SLwork) { }
   ~SmoLenFindWorker() { }*/

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
            const pnodPtrT pnod = static_cast<pnodPtrT>(curPart);
            const _partT* partPtr = static_cast<_partT*>(pnod->partPtr);

            const cType noneigh = partPtr->noneigh;
            const fType smass  = static_cast<fType>(noneigh);
            const fType srad   = BHTreeWorker::maxMassEncloseRad(pnod, smass);
            /*std::list<_partT*> neighs =
               NeighFindWorker<_partT>::operator()(pnod, srad);*/

            parentT::Func.neighList.clear();
            parentT::neighExecFunc(pnod, srad);
            //std::list<PartPtrDistT>& neighList(parentT::Func.neighList);

      //return(parentT::Func.neighList);


            //NeighDistComp myNeighDistComp(partPtr);
            //neighs.sort(myNeighDistComp);

            //if ( neighs.size() <= noneigh )

         }
         curPart = curPart->next;
      }
   }

private:
   // true if _p1 is nearer to _ref than to _p2
   class NeighDistComp
   {
public:
      NeighDistComp(_partT* _ref) : ref(_ref->pos) { }
      ~NeighDistComp() { }

      bool operator()(_partT* _p1, _partT* _p2)
      {
         const vect3dT rvec1 = ref - _p1->pos;
         const vect3dT rvec2 = ref - _p2->pos;

         return(dot(rvec1, rvec1) < dot(rvec2, rvec2));
      }

private:
      const vect3dT ref;
   };
};
};


#endif

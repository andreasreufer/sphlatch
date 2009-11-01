#ifndef BHTREE_WORKER_GRAV_CPP
#define BHTREE_WORKER_GRAV_CPP

/*
 *  bhtree_worker_grav.cpp
 *
 *  Created by Andreas Reufer on 14.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"
#include "bhtree.h"

namespace sphlatch {
template<typename _MAC, typename _partT>
class GravityWorker : public BHTreeWorkerRO {
public:
   GravityWorker(const treePtrT _treePtr) : BHTreeWorkerRO(_treePtr) { }
   GravityWorker(const GravityWorker& _gw) : BHTreeWorkerRO(_gw) { }
   ~GravityWorker() { }

   void calcGravity(const czllPtrT _czll);

private:
   _MAC mac;

   void calcGravPart(const pnodPtrT _part);

   void interactPartCell();
   void interactPartPart();

   vect3dT acc, ppos;
};

template<typename _MAC, typename _partT>
void GravityWorker<_MAC, _partT>::calcGravity(const czllPtrT _czll)
{
  nodePtrT curPart = _czll->chldFrst;
  const nodePtrT stopChld = _czll->chldLast->next;

  while ( curPart != stopChld )
  {
    calcGravPart(curPart);
    curPart = curPart->next;
  }
};

template<typename _MAC, typename _partT>
void GravityWorker<_MAC, _partT>::calcGravPart(const pnodPtrT _part)
{
   nodePtrT const curPartPtr = _part;

   ppos = _part->pos;
   acc  = 0, 0, 0;

   ///
   /// the complete tree walk
   ///
   goRoot();
   do
   {
      if (not curPtr->isParticle)
      {
         if (mac(static_cast<gcllPtrT>(curPtr), curPartPtr))
         {
            interactPartCell();
            goSkip();
         }
         else
            goNext();
      }
      else
      {
         if (curPtr != curPartPtr)
         {
            interactPartPart();
         }
         goNext();
      }
   } while (curPtr != NULL);

   static_cast<_partT*>(_part->partPtr)->acc += acc;
}

template<typename _MAC, typename _partT>
void GravityWorker<_MAC, _partT>::interactPartPart()
{
   //FIXME: try using blitz++ functions
   const fType rx = ppos[0] - static_cast<pnodPtrT>(curPtr)->pos[0];
   const fType ry = ppos[1] - static_cast<pnodPtrT>(curPtr)->pos[1];
   const fType rz = ppos[2] - static_cast<pnodPtrT>(curPtr)->pos[2];
   const fType rr = rx * rx + ry * ry + rz * rz;
   const fType r  = sqrt(rr);

   const fType m = static_cast<pnodPtrT>(curPtr)->m;

   const fType mOr3 = m / (rr * r);

   acc[0] -= mOr3 * rx;
   acc[1] -= mOr3 * ry;
   acc[2] -= mOr3 * rz;
}

template<typename _MAC, typename _partT>
void GravityWorker<_MAC, _partT>::interactPartCell()
{
   const fType rx = mac.rx;
   const fType ry = mac.ry;
   const fType rz = mac.rz;
   const fType rr = mac.rr;
   const fType r  = sqrt(rr);

   const fType Or3 = 1. / (r * rr);

   const fType m = static_cast<qcllPtrT>(curPtr)->m;

   acc[0] -= m * Or3 * rx;
   acc[1] -= m * Or3 * ry;
   acc[2] -= m * Or3 * rz;

   const fType Or5 = Or3 / rr;
   const fType Or7 = Or5 / rr;

   const fType q11 = static_cast<qcllPtrT>(curPtr)->q11;
   const fType q22 = static_cast<qcllPtrT>(curPtr)->q22;
   const fType q33 = static_cast<qcllPtrT>(curPtr)->q33;
   const fType q12 = static_cast<qcllPtrT>(curPtr)->q12;
   const fType q13 = static_cast<qcllPtrT>(curPtr)->q13;
   const fType q23 = static_cast<qcllPtrT>(curPtr)->q23;

   const fType q1jrj   = q11 * rx + q12 * ry + q13 * rz;
   const fType q2jrj   = q12 * rx + q22 * ry + q23 * rz;
   const fType q3jrj   = q13 * rx + q23 * ry + q33 * rz;
   const fType qijrirj = q11 * rx * rx +
                         q22 * ry * ry +
                         q33 * rz * rz +
                         2. * q12 * rx * ry +
                         2. * q13 * rx * rz +
                         2. * q23 * ry * rz;

   acc[0] += (Or5) * (q1jrj) - (Or7) * (2.5 * qijrirj * rx);
   acc[1] += (Or5) * (q2jrj) - (Or7) * (2.5 * qijrirj * ry);
   acc[2] += (Or5) * (q3jrj) - (Or7) * (2.5 * qijrirj * rz);
}

class fixThetaMAC {
   typedef monopoleCellNode*     mcllPtrT;
   typedef quadrupoleCellNode*   qcllPtrT;
   typedef particleNode*         pnodPtrT;
public:
   fType rx, ry, rz, rr;
   bool operator()(const qcllPtrT _cell, const pnodPtrT _part)
   {
      const fType theta  = 0.6;
      const fType clsz   = _cell->clSz;
      const fType theta2 = theta * theta;

      rx = _cell->com[0] - _part->pos[0];
      ry = _cell->com[1] - _part->pos[1];
      rz = _cell->com[2] - _part->pos[2];
      rr = rx * rx + ry * ry + rz * rz;

      return((clsz * clsz / rr) < theta2);
   }
};
};

#endif

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
#include "timer.cpp"

namespace sphlatch {
template<typename _macT, typename _partT>
/*class GravityWorker : public BHTreeWorker {
public:
   GravityWorker(const treePtrT _treePtr,
                 const fType    _G) : BHTreeWorker(_treePtr), G(_G) { }
   GravityWorker(const GravityWorker& _gw) : BHTreeWorker(_gw), G(_gw.G) { }
   ~GravityWorker() { }*/

class GravityWorker : public BHTreeWorkerRO {
public:
   GravityWorker(const treePtrT _treePtr,
                 const fType    _G) : BHTreeWorkerRO(_treePtr), G(_G) { }
   GravityWorker(const GravityWorker& _gw) : BHTreeWorkerRO(_gw), G(_gw.G) { }
   ~GravityWorker() { }

   void calcGravity(const czllPtrT _czll);
   void calcGravPart(const pnodPtrT _part);

   typedef sphlatch::Timer   timerT;

private:
   _macT  MAC;
   timerT Timer;

   void calcGravPartAlt(const pnodPtrT _part);
   void calcGravRec();

   void interactPartCell();
   void interactPartPart();

   vect3dT  acc, ppos;
   nodePtrT recCurPartPtr;

   fType splineOSmoR3(const fType _r, const fType _h);

protected:
   const fType G;
};

template<typename _macT, typename _partT>
void GravityWorker<_macT, _partT>::calcGravity(const czllPtrT _czll)
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
         calcGravPart(static_cast<pnodPtrT>(curPart));
      curPart = curPart->next;
   }
   const double compTime = Timer.getRoundTime();
   _czll->compTime += static_cast<fType>(compTime);
}

template<typename _macT, typename _partT>
void GravityWorker<_macT, _partT>::calcGravPart(const pnodPtrT _part)
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
         if (MAC(static_cast<qcllPtrT>(curPtr),
                 static_cast<pnodPtrT>(curPartPtr)))
         {
            interactPartCell();
            goSkip();
         }
         else
         {
            goNext();
         }
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

   static_cast<_partT*>(_part->partPtr)->acc += G * acc;
}

template<typename _macT, typename _partT>
void GravityWorker<_macT, _partT>::calcGravPartAlt(const pnodPtrT _part)
{
   ppos          = _part->pos;
   acc           = 0, 0, 0;
   recCurPartPtr = _part;
   ///
   /// the complete tree walk
   ///
   goRoot();
   calcGravRec();

   static_cast<_partT*>(_part->partPtr)->acc += G * acc;
}

template<typename _macT, typename _partT>
void GravityWorker<_macT, _partT>::calcGravRec()
{
   if (curPtr->isParticle)
   {
      if (curPtr != recCurPartPtr)
         interactPartPart();
   }
   else
   {
      if (MAC(static_cast<qcllPtrT>(curPtr),
              static_cast<pnodPtrT>(recCurPartPtr)))
      {
         interactPartCell();
      }
      else
      {
         for (size_t i = 0; i < 8; i++)
         {
            if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
            {
               goChild(i);
               calcGravRec();
               goUp();
            }
         }
      }
   }
}

template<typename _macT, typename _partT>
void GravityWorker<_macT, _partT>::interactPartPart()
{
   //FIXME: try using blitz++ functions
   const fType rx = ppos[0] - static_cast<pnodPtrT>(curPtr)->pos[0];
   const fType ry = ppos[1] - static_cast<pnodPtrT>(curPtr)->pos[1];
   const fType rz = ppos[2] - static_cast<pnodPtrT>(curPtr)->pos[2];
   const fType m  = static_cast<pnodPtrT>(curPtr)->m;
   const fType rr = rx * rx + ry * ry + rz * rz;
   const fType r  = sqrt(rr);

#ifdef SPHLATCH_GRAVITY_SPLINESMOOTHING
   const fType h =
      static_cast<_partT*>(static_cast<pnodPtrT>(curPtr)->partPtr)->h;
   const fType mOr3 = m * splineOSmoR3(r, h);
#elif SPHLATCH_GRAVITY_EPSSMOOTHING
   const fType eps =
      static_cast<_partT*>(static_cast<pnodPtrT>(curPtr)->partPtr)->eps;
   const fType re   = r + eps;
   const fType mOr3 = m / (re * re * re);
#else
   const fType mOr3 = m / (rr * r);
#endif

   acc[0] -= mOr3 * rx;
   acc[1] -= mOr3 * ry;
   acc[2] -= mOr3 * rz;
}

template<typename _macT, typename _partT>
void GravityWorker<_macT, _partT>::interactPartCell()
{
   //FIXME: check if fetching those values again is less costly
   const fType rx = MAC.rx;
   const fType ry = MAC.ry;
   const fType rz = MAC.rz;

   const fType rr = MAC.rr;
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

///
/// gravitational spline softening for B Spline kernel from
/// Hernquist & Katz 1989
///
template<typename _macT, typename _partT>
fType GravityWorker<_macT, _partT>::splineOSmoR3(const fType _r, const fType _h)
{
   const fType u = _r / _h;

   if (u >= 2.)
   {
      const fType r3 = _r * _r * _r;
      return(1. / r3);
   }
   else
   if (u > 1.)
   {
      const fType r3 = _r * _r * _r;
      return((1. / r3) * (
                -(1. / 15.)
                + (8. / 3.) * u * u * u
                - 3. * u * u * u * u
                + (6. / 5.) * u * u * u * u * u
                - (1. / 6.) * u * u * u * u * u * u
                ));
   }
   else
   {
      const fType h3 = _h * _h * _h;
      return((1. / h3) * (
                +(4. / 3.)
                - (6. / 5.) * u * u
                + (1. / 2.) * u * u * u
                ));
   }
}

class fixThetaMAC {
   typedef monopoleCellNode*     mcllPtrT;
   typedef quadrupoleCellNode*   qcllPtrT;
   typedef particleNode*         pnodPtrT;
public:
   fType rx, ry, rz, rr;
   bool operator()(const qcllPtrT _cell, const pnodPtrT _part)
   {
      const fType theta  = 0.70;
      const fType clsz   = _cell->clSz;
      const fType theta2 = theta * theta;

      rx = _part->pos[0] - _cell->com[0];
      ry = _part->pos[1] - _cell->com[1];
      rz = _part->pos[2] - _cell->com[2];
      rr = rx * rx + ry * ry + rz * rz;

      return((clsz * clsz / rr) < theta2);
   }
};
};

#endif

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
#include "bhtree_particle.h"
#include "bhtree.h"

// FIXME: use function template for MAC

namespace sphlatch {
template<typename _MAC>
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



};

template<typename _MAC>
void GravityWorker<_MAC>::calcGravPart(const pnodPtrT _part)
{
   nodePtrT const curPartPtr = _part;

   fType accX = 0., accY = 0., accZ = 0.;

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
            //calcGravMP()
            goSkip();
         }
         else
            goNext();
      }
      else
      {
         if (curPtr != curPartPtr)
          //calcGravPart()
         goNext();
      }
   } while (curPtr != NULL);
}


class fixThetaMAC {
   typedef monopoleCellNode*     mcllPtrT;
   typedef quadrupoleCellNode*   qcllPtrT;
   typedef particleNode*         pnodPtrT;
public:
   bool operator()(const qcllPtrT _cell, const pnodPtrT _part)
   {
      const fType theta = 0.6;

      const fType rx     = _cell->com[0] - _part->pos[0];
      const fType ry     = _cell->com[1] - _part->pos[1];
      const fType rz     = _cell->com[2] - _part->pos[2];
      const fType rr     = rx * rx + ry * ry + rz * rz;
      const fType clsz   = _cell->clSz;
      const fType theta2 = theta*theta;

      return((clsz * clsz / rr) < theta2);
   }
};
};

#endif

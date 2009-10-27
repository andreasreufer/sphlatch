#ifndef BHTREE_WORKER_GRAV_CPP
#define BHTREE_WORKER_GRAV_CPP

/*
 *  bhtree_worker_grav.cpp
 *
 *  Created by Andreas Reufer on 14.12.08.
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.h"

// FIXME: use function template for MAC

namespace sphlatch {
class BHTreeWorkerGrav : public BHTreeWorkerRO {
public:
   BHTreeWorkerGrav(treePtrT _treePtr) : BHTreeWorkerRO(_treePtr)
   {
      //const size_t tid = omp_get_thread_num();
      //std::cout << tid << ":tpc:" << curPtr << ":" << this << "\n";
   }

   BHTreeWorkerGrav(const BHTreeWorkerGrav& _gravwork)
      : BHTreeWorkerRO(_gravwork)
   {
      //const size_t tid = omp_get_thread_num();
      //std::cout << tid << ":cpc:" << curPtr << ":" << this << "\n";
   }

   ~BHTreeWorkerGrav() { }

private:
   void calcGravParticle(const size_t _i);
};


void BHTreeWorkerGrav::calcGravParticle(const size_t _i)
{
   /*const fType posX = pos(_i, X);
   const fType posY = pos(_i, Y);
   const fType posZ = pos(_i, Z);*/

   //nodePtrT const curPartPtr = treePtr->partProxies[_i];
   nodePtrT const curPartPtr = NULL;

   fType accX = 0., accY = 0., accZ = 0.;

   ///
   /// the complete tree walk
   ///
   goRoot();
   do
   {
      if (not curPtr->isParticle)
      {
         if (true) // MAC()
         {
            //calcGravMP
            goSkip();
         }
         else
            goNext();
      }
      else
      {
         if (curPtr != curPartPtr)
         {
         //calcGravPart
         }
         goNext();
      }
   } while (curPtr != NULL);

   /*acc(_i, X) += accX;
   acc(_i, Y) += accY;
   acc(_i, Z) += accZ;*/
}
};

#endif

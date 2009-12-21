#ifndef BHTREE_WORKER_MP_CPP
#define BHTREE_WORKER_MP_CPP

/*
 *  bhtree_worker_mp.cpp
 *
 *  Created by Andreas Reufer on 27.10.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"

namespace sphlatch {
class BHTreeMPWorker : public BHTreeWorker {
public:

   BHTreeMPWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeMPWorker(const BHTreeMPWorker& _mp) : BHTreeWorker(_mp) { }
   ~BHTreeMPWorker() { }

   void calcMultipoles(const czllPtrT _czllPtr);
   void calcMultipolesCZ();

private:

   void MPRec();
   void MPRecCZ();
};

void BHTreeMPWorker::calcMultipoles(const czllPtrT _czllPtr)
{
   curPtr = _czllPtr;
   MPRec();
}

void BHTreeMPWorker::MPRec()
{
   if (not curPtr->isParticle)
   {
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            if (not static_cast<gcllPtrT>(curPtr)->child[i]->isParticle)
            {
               goChild(i);
               MPRec();
               goUp();
            }
         }
      }
      static_cast<qcllPtrT>(curPtr)->calcMultipole();
   }
}

void BHTreeMPWorker::calcMultipolesCZ()
{
   goRoot();
   MPRecCZ();
}

void BHTreeMPWorker::MPRecCZ()
{
   if (not curPtr->isParticle)
   {
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
         // omit the particle check, as there should be no particles
         // in the CZ part of the tree
            if (not static_cast<gcllPtrT>(curPtr)->child[i]->atBottom)
            {
               goChild(i);
               MPRecCZ();
               goUp();
            }
         }
      }
      static_cast<qcllPtrT>(curPtr)->calcMultipole();
   }
}
};

#endif

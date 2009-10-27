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
   BHTreeMPWorker(const BHTree& _mp) : BHTreeWorker(_mp) { }
   ~BHTreeMPWorker() { }

   void calcMultipoles(const czllPtrT _czllPtr);
   void calcMultipolesCZ();

private:

   void MPRec();
   void MPRecCZ();
};

void BHTreeMPWorker::calcMultipoles(const czllPtrT _czllPtr)
{
  
}

void MPRec()
{ }

void BHTreeMPWorker::calcMultipolesCZ()
{ }

void MPRecCZ()
{ }
};

#endif

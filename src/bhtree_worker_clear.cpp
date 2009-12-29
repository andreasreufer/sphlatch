#ifndef BHTREE_WORKER_CLEAR_CPP
#define BHTREE_WORKER_CLEAR_CPP

/*
 *  bhtree_worker_clear.cpp
 *
 *  Created by Andreas Reufer on 29.12.09
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"

namespace sphlatch {
class BHTreeClear : public BHTreeWorker {
public:
   BHTreeClear(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   ~BHTreeClear() { }

   void operator()();

private:
   void recursor();
};

void BHTreeClear::operator()()
{
  goRoot();
};

void BHTreeClear::recursor()
{
   if (not curPtr->isParticle)
   {
      for (size_t i = 0; i < 8; i++)
      {
         goChild(i);
         recursor();
         goRoot();
      }
      if (curPtr->isCZ)
         delete static_cast<czllPtrT>(curPtr);
      else
         delete static_cast<qcllPtrT>(curPtr);
   }
   else
      delete static_cast<pnodPtrT>(curPtr);
};


};
#endif

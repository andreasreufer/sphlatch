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

   void operator()(const nodePtrT _node);

private:
   void recursor();
};

void BHTreeClear::operator()(const nodePtrT _node)
{
   curPtr = _node;
   recursor();
}

void BHTreeClear::recursor()
{
   if (not curPtr->isParticle)
   {
      for (size_t i = 0; i < 8; i++)
      {
         if (static_cast<gcllPtrT>(curPtr)->child[i] != NULL)
         {
            if (static_cast<gcllPtrT>(curPtr)->child[i]->isParticle)
               delete static_cast<pnodPtrT>(
                  static_cast<gcllPtrT>(curPtr)->child[i]);
            else
            {
               goChild(i);
               recursor();
               goUp();
               if (static_cast<gcllPtrT>(curPtr)->child[i]->isCZ)
                  delete static_cast<czllPtrT>(
                     static_cast<gcllPtrT>(curPtr)->child[i]);
               else
                  delete static_cast<qcllPtrT>(
                     static_cast<gcllPtrT>(curPtr)->child[i]);
            }
            static_cast<gcllPtrT>(curPtr)->child[i] = NULL;
         }
      }
   }
}
};
#endif

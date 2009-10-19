#ifndef BHTREE_HOUSEKEEPER_H  
#define BHTREE_HOUSEKEEPER_H  

/*
 *  bhtree_housekeeper.h
 *
 *  Created by Andreas Reufer on 19.10.09
 *  Copyright 2007 University of Berne. All rights reserved.
 *
 */

#include <vector>

#include "bhtree_worker.h"

namespace sphlatch {
class BHTreeHousekeeper : public BHTreeWorker {
public:

   BHTreeHousekeeper(const treePtrT _treePtr);
   BHTreeHousekeeper(const BHTreeHousekeeper& _hk);
   ~BHTreeHousekeeper();

   ///
   /// set the next & skip pointers
   ///
   void setNext(const czllPtrT _czll);
   void setNextCZ();
   void setSkip();

private:
   void setNextRecursor();
   void setNextCZRecursor();

   nodePtrT lastPtr;

   std::vector<gcllPtrT> lastSkipeeAtDepth;
};

};

#endif

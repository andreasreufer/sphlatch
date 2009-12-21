#ifndef BHTREE_WORKER__CPP
#define BHTREE_WORKER__CPP

/*
 *  bhtree_worker_.cpp
 *
 *  Created by Andreas Reufer on 14.12.08
 *  Copyright 2009 University of Berne. All rights reserved.
 *
 */

#include "bhtree_worker.cpp"

namespace sphlatch {
class BHTreeMyWorker : public BHTreeWorker {
public:

   BHTreeMyWorker(const treePtrT _treePtr) : BHTreeWorker(_treePtr) { }
   BHTreeMyWorker(const BHTree& _) : BHTreeWorker(_) { }
   ~BHTreeMyWorker() { }
};
};

#endif
